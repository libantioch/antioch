//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
#include "antioch/ReactionChem.hpp"
#include <algorithm>
#include <iomanip>
#include <sstream>

namespace Antioch{
ReactionChem::ReactionChem(const ReactionChem &rhs)
{
  setGlobalRateConstant(rhs.getGlobalRate());
  setAllReactiveChannels(rhs.getAllReactiveChannels());
  id = rhs.getID();
}

void ReactionChem::setGlobalRateConstant(const KineticsChem &chem)
{
   globalRate.setNonReactivePart(chem.getSimpleChemPart());
   globalRate.setKineticsOnly(chem);
}

void ReactionChem::setAllReactiveChannels(const std::vector<SimpleChem> &chem)
{
//equalizing nchannels if needed
  for(int iout = getNChannels(); iout >= (int)chem.size() && iout > 0; iout--)
  {
    reactiveChannel.pop_back();
  }

  for(unsigned int ichan = 0; ichan < chem.size(); ichan++)
  {
    if((int)ichan < getNChannels())
    {
      setReactiveChannel(chem[ichan],ichan);
    }else
    {
      addReactiveChannel(chem[ichan]);
    }
  }

  branchingRatios.resize(chem.size());
  makeBranchingRatios();
}

void ReactionChem::setReaction(const ReactionChem &rhs)
{
  setGlobalRateConstant(rhs.getGlobalRate());
  setAllReactiveChannels(rhs.getAllReactiveChannels());
  id = rhs.getID();
}

void ReactionChem::setKinetics(const std::string &kinProc, const std::string &kinModel, const std::vector<ParameterPhy> &pars)
{
  globalRate.setKineticsProcess(kinProc);
  globalRate.setKineticsModel(kinModel);
  globalRate.storeParameters(pars);
  globalRate.setReactionFromStorage();
}

bool ReactionChem::GotTheseProducts(const std::vector<std::string> &mols)  const
{
  bool out(false);
  for(unsigned int i = 0; i < reactiveChannel.size(); i++)
  {
    out = out || reactiveChannel[i].areSameMolecules(mols);
  }
  return out;
}

int ReactionChem::WhereTheseProducts(const std::vector<std::string> &mols)  const
{
  int out(-1);
  for(unsigned int i = 0; i < reactiveChannel.size(); i++)
  {
    if(reactiveChannel[i].areSameMolecules(mols))out = (int)i;
  }
  return out;
}

int ReactionChem::nPars() const
{
  int n = globalRate.getNPar();
  for(unsigned int i = 0 ; i < reactiveChannel.size(); i++)
  {
    n += reactiveChannel[i].getNPar();
  }
  return n;  
}

void ReactionChem::densifyParametersToThis(const ParameterPhy &par)  
{
  for(int ipg = 0; ipg < globalRate.getNPar(); ipg++)
  {
    globalRate.densifyParametersToThis(par);
  }
  for(unsigned int i = 0 ; i < reactiveChannel.size(); i++)
  {
    reactiveChannel[i].densifyParametersToThis(par);
  }
}

void ReactionChem::densifyAllUnitsIntra()
{
  globalRate.densifyParametersUnit();
  for(unsigned int i = 0 ; i < reactiveChannel.size(); i++)
  {
    reactiveChannel[i].densifyParametersUnit();
  }
}

void ReactionChem::densifyAllUnitsExtra()
{
  globalRate.densifyParametersUnit();
  for(unsigned int i = 0 ; i < reactiveChannel.size(); i++)
  {
    int k = i - 1;
    while(k >= -1)
    {
      if(k < 0)
      {
        reactiveChannel[i].densifyParametersUnitExternal(globalRate);
      }else
      {
        reactiveChannel[i].densifyParametersUnitExternal(reactiveChannel[k]);
      }
      k--;
    }
  }
}

void ReactionChem::densifyAllUnitsExternal(ReactionChem &out)  
{
  globalRate.densifyParametersUnitExternal(out.getGlobalRate());
  for(int i = 0 ; i < (int)reactiveChannel.size(); i++)
  {
/* k == 1 : out.globalRate
 * 1 < k <= out.getNChannels() + 1      : out.branchingRatios
 * k > out.getNChannels() + 1 && i == 0 : this->globalRate
 * k > out.getNChannels() + 1 && i > 0  : this->branchingRatios
 */
    int k = (out.getNChannels() + 1) + (i + 2); //(out channels + out global rate) + (this channels + this global rate)
    while(k > 0)
    {
      if(k > out.getNChannels() + 1)//intra ReactionChem
      {
        if(i == 0)
        {
          reactiveChannel[i].densifyParametersUnitExternal(globalRate);
        }else
        {
          reactiveChannel[i].densifyParametersUnitExternal(reactiveChannel[i - 1]);
        }
      }else
      {// extra ReactionChem
        if(k == 1)
        {
          reactiveChannel[i].densifyParametersUnitExternal(out.getGlobalRate());
        }else
        {
          reactiveChannel[i].densifyParametersUnitExternal(out.getReactiveChannel(k - 2));
        }
      }
      k--;
    }
  }
}


//default branching ratio => br = 1, dbr = 0 (none)
void ReactionChem::addAbsoluteOneBranchingRatio(const std::vector<std::string> &mols)
{
  SimpleChem tmpBr;
  tmpBr.addParameter(ParameterPhy("branching ratio",1.,0.,CORE_UNCERTAINTY_TYPE_NONE,std::string()));
  tmpBr.setDescription("constant model: par = constant");
  tmpBr.setMolecules(mols);
  reactiveChannel.push_back(tmpBr);
}

void ReactionChem::ReactiveChannelsClearAll()
{
  for(int ichan = 0; ichan < getNChannels(); ichan++)
  {
    reactiveChannel[ichan].clear();
  }
}

void ReactionChem::makeBranchingRatio(int ichan)
{
  if((int)branchingRatios.size() <= ichan)
  {
      branchingRatios.resize(ichan + 1);
  }

  branchingRatios[ichan] = ParameterPhy("branching ratio",1.0,0.,CORE_UNCERTAINTY_TYPE_ABSOLUTE,"");

  for(int ilevel = 0; ilevel < ReactiveChannelNParameters(ichan); ilevel++)
  {
      branchingRatios[ichan] *= getReactiveChannelParameter(ilevel,ichan);
  }
}

void ReactionChem::makeBranchingRatios()
{
  for(int ichan = 0; ichan < (int)reactiveChannel.size(); ichan++)
  {
     makeBranchingRatio(ichan);
  }

}

void ReactionChem::ReactiveChannelsClearParameters()
{
  for(int ichan = 0; ichan < getNChannels(); ichan++)
  {
    ReactiveChannelClearParameters(ichan);
  }
}

void ReactionChem::ReactiveChannelsClearParametersValues()
{
  for(int ichan = 0; ichan < getNChannels(); ichan++)
  {
    ReactiveChannelClearParametersValues(ichan);
  }
}

void ReactionChem::clearAll()
{
  GlobalRateConstantClear();
  ReactiveChannelsEraseAll();
}

void ReactionChem::resizeAll(int newSize)
{
  for(int ip = 0; ip < GlobalRateConstantNParameters(); ip++)
  {
    GlobalRateConstantResize(newSize,ip);
  }
  for(int ichan = 0; ichan < getNChannels(); ichan++)
  {
    for(int ipar = 0; ipar < reactiveChannel[ichan].getNPar(); ipar++)
    {
      ReactiveChannelResize(newSize,ipar,ichan);
    }
  }
}

const std::string ReactionChem::report(int ichan)
{
  std::ostringstream outstream;

  outstream << std::endl;
  for(int i =0; i < NReactants() - 1; i++)
  {
    outstream << Reactant(i) << " + ";
  }
  outstream << Reactant(NReactants() - 1) << " -> ";
  for(int i =0; i < NProducts(ichan) - 1; i++)
  {
    outstream << Product(i,ichan) << " + ";
  }
  outstream << Product(NProducts(ichan) - 1,ichan) << std::endl;

  outstream << getGlobalRateConstantKineticsProcess() << ": " << getGlobalRateConstantKineticsModel() << " kinetics model" << std::endl;

  outstream << std::setw(32) << std::left << "parameter"
            << "\n\t"
            << std::setw(30) << "unit"
            << std::setw(15) << "q_05"
            << std::setw(15) << "q_95"
            << std::setw(15) << "q_50"
            << std::setw(15) << "mean"
            << std::setw(15) << "std dev"
            << std::endl;

  for(int ipar = 0; ipar < GlobalRateConstantNParameters(); ipar++)
  {
    std::vector<double> par = GlobalRateConstantParValues(ipar);
    if(GlobalRateConstantParName(ipar).find("PreExp") != std::string::npos)
    {
        for(int ival = 0; ival < GlobalRateConstantParNValues(ipar); ival++)
        {
            par[ival] *= ReactiveChannelBranchingRatioValue(ival,ichan);
        }
    }

// quantiles
    sort(par.begin(),par.end());
    double q05,q95,q50;
    if(par.size() % 20 == 0)
    {
       q05 = par[par.size()/20];
       q95 = par[19*par.size()/20];
    }else
    {
       q05 = (par[par.size()/20] + par[par.size()/20 + 1] ) / 2.;
       q95 = (par[19*par.size()/20] + par[19*par.size()/20 + 1] ) / 2.;
    }
    (par.size() % 2 == 0)?
                q50 = par[par.size()/2]:
                q50 = par[par.size()/2 + 1];
//mean & var
    double mean(0.),var(0.);
    for(unsigned int i = 0; i < par.size(); i++)
    {
      mean += par[i];
      if(i > 0)
      {
        double vf = (double)i + 1.;
        double rnnf = 1./(vf * (vf - 1.));
        var += rnnf * (vf * par[i] - mean) * (vf * par[i] - mean);
      }
    }
    mean /= (double)par.size();
    var /= (double)par.size();

    outstream << std::setw(30) << std::left << GlobalRateConstantParName(ipar) 
              << "\n\t"   
              << std::setw(30) << getGlobalRateConstantParUnit(ipar)
              << std::setw(15) << q05
              << std::setw(15) << q95
              << std::setw(15) << q50
              << std::setw(15) << mean
              << std::setw(15) << std::sqrt(var)
              << std::endl;
  }

  return outstream.str();
  
}

const std::string ReactionChem::reportAll()
{
  std::string out("");
  for(int i = 0; i < getNChannels(); i++)
  {
    out += report(i) + "\n";
  }

  return out;
}

const std::vector<ParameterPhy> ReactionChem::ReactiveChannelBranchingRatios()
{
  makeBranchingRatios();

  return branchingRatios;
}

ParameterPhy ReactionChem::ReactiveChannelBranchingRatio(int ichan)
{
  makeBranchingRatio(ichan);

  return branchingRatios[ichan];
}

double ReactionChem::ReactiveChannelBranchingRatioValue(int ival, int ichan)
{
//br model (if it exists) not yet supported
  makeBranchingRatio(ichan);

  return branchingRatios[ichan].value(ival);
}

const std::string ReactionChem::ReactionEquation() const
{
  std::string equa;
  for(int ireac = 0; ireac < NReactants(); ireac++)
  {
    equa += Reactant(ireac);
    (ireac == NReactants() - 1)?equa += " ->":equa += " + ";
  }
  equa += "\n";
  for(int ichan = 0; ichan < getNChannels(); ichan++)
  {
    equa += "\t";
    for(int iprod = 0; iprod < NProducts(ichan); iprod++)
    {
       equa += Product(iprod,ichan);
      (iprod == NProducts(ichan) - 1)?equa += "\n":equa += " + ";
    }
  }

  return equa;
}

const std::string ReactionChem::ReactionChannelEquation(int ichan) const
{
  std::string equa;
  for(int ireac = 0; ireac < NReactants(); ireac++)
  {
    equa += Reactant(ireac);
    (ireac == NReactants() - 1)?equa += " -> ":equa += " + ";
  }
  for(int iprod = 0; iprod < NProducts(ichan); iprod++)
  {
     equa += Product(iprod,ichan);
    (iprod == NProducts(ichan) - 1)?equa += " ":equa += " + ";
  }

  return equa;

}

void ReactionChem::GlobalRateConstantShowAll(std::ostream &out) const
{
  out << "Non reactive part" << std::endl;
  globalRate.showAll(out);
  out << "Reactive part" << std::endl;
  globalRate.showKinetics(out);
}

// showAll()
const void ReactionChem::showAll(std::ostream &out) const
{

  for(int i = 0; i < 100; i++)out << "-";
  out << "\nDescription of a ReactionChem object." << std::endl;
  out << "\nThis reaction has " << getNChannels() << " reactive channels" << std::endl;

  ShowReactionResume(out);

  GlobalRateConstantShowAll(out);
  for(int i = 0; i < getNChannels(); i++)
  {
    ReactiveChannelShowAll(i,out);
  }

  out << "End of the description of the ReactionChem object." << std::endl;
  for(int i = 0; i < 100; i++)out << "-";
  out << std::endl;

}
}
