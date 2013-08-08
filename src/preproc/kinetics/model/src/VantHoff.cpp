#include "antioch/VantHoff.hpp"
#include "antioch/phymath.hpp"

namespace Antioch{

VantHoff::VantHoff(const std::vector<ParameterPhy> &pars)
{
  init(pars);
}

void VantHoff::init(const std::vector<ParameterPhy> &pars)
{
  std::vector<const ParameterPhy*> spars = whatParametersAreWe(pars);

  setEa(*spars[0]);
  setD(*spars[1]);
  setBeta(*spars[2]);
  setT0(*spars[3]);
  setPreExp(*spars[4]);
}

void VantHoff::reset(const std::vector<double> &pars)
{
  const std::string method("void Kooij::reset(const std::vector<double> &)");

//number of parameters
  if(pars.size() != nPars())antiochError(method,"Not the right number of parameters");
  resetPreExp(pars[0]);
  resetBeta(pars[1]);
  resetEa(pars[2]);
  resetD(pars[3]);
}

void VantHoff::checkParams() const
{
  checkA("VantHoff::checkParams() const");
  checkParam(beta,"beta","VantHoff::checkParams() const");
  checkParam(Ea,"Ea","VantHoff::checkParams() const");
  checkParam(D,"D","VantHoff::checkParams() const");
}

const ParameterPhy VantHoff::rateConstant(int nk) 
{
  checkT("VantHoff");
  checkParams();

  resetPar(ACal,A.value(nk));
  resetPar(betaCal,beta.value(nk));
  resetPar(expPaCal,expPa.value(nk));
  resetPar(DCal,D.value(nk));

  ParameterPhy out = kEquation<ParameterPhy>(*T,ACal,betaCal,expPaCal,DCal,T0);
  out.setNamePar("VantHoff rate constant");
  if(T->isUncDef()) //T is uncertain
        out.correcUncAddBiais(2. * pow(out,2) * (DCal*betaCal/(*T) + 
                                                 DCal*expPaCal/(pow(*T,2)) + 
                                                 betaCal*expPaCal/(pow(*T,3))),*T);
  return out;
}

double VantHoff::getRateConstantT(double Tem, int i) const
{
  return kEquation<double>(Tem,A.value(i),beta.value(i), expPa.value(i), D.value(i), T0.value());
}

const std::vector<ParameterPhy> VantHoff::getParameters() const
{
  std::vector<ParameterPhy> out;
  out.push_back(A);
  out.push_back(beta);
  out.push_back(Ea);
  out.push_back(D);
  out.push_back(T0);

  return out;
}

const std::string VantHoff::getDefaultUnit(int i) const
{
  const std::string method("const std::string VantHoff::getDefaultUnit(int ) const");
  switch(i)
  {
     case 0:
     {
       return getDefaultPreExpUnit();
       break;
     }
     case 1:
     {
       return getDefaultBetaUnit();
       break;
     }
     case 2:
     {
       return getDefaultEaUnit();
       break;
     }
     case 3:
     {
       return getDefaultDUnit();
       break;
     }
  }
  antiochError(method,"Asking for a non existing parameter");

  return std::string();
}

double VantHoff::getDefaultMin(int i) const
{
  const std::string method("double VantHoff::getDefaultMin(int ) const");
  switch(i)
  {
     case 0:
     {
       return getDefaultPreExpMin();
       break;
     }
     case 1:
     {
       return getDefaultBetaMin();
       break;
     }
     case 2:
     {
       return getDefaultEaMin();
       break;
     }
     case 3:
     {
       return getDefaultDMin();
       break;
     }
  }
  antiochError(method,"Asking for a non existing parameter");

  return 0.;
}

double VantHoff::getDefaultMax(int i) const
{
  const std::string method("double VantHoff::getDefaultMax(int ) const");
  switch(i)
  {
     case 0:
     {
       return getDefaultPreExpMax();
       break;
     }
     case 1:
     {
       return getDefaultBetaMax();
       break;
     }
     case 2:
     {
       return getDefaultEaMax();
       break;
     }
     case 3:
     {
       return getDefaultDMax();
       break;
     }
  }
  antiochError(method,"Asking for a non existing parameter");

  return 0.;
}

const ParameterPhy VantHoff::getDefaultParameter(int i) const
{
  const std::string method("const ParameterPhy VantHoff::getDefaultParameter(int )");
  switch(i)
  {
     case 0:
     {
       return getDefaultPreExp();
       break;
     }
     case 1:
     {
       return getDefaultBeta();
       break;
     }
     case 2:
     {
       return getDefaultEa();
       break;
     }
     case 3:
     {
       return getDefaultD();
       break;
     }
  }
  antiochError(method,"Asking for a non existing parameter");

  return ParameterPhy();
}

std::vector<const ParameterPhy*> VantHoff::whatParametersAreWe(const std::vector<ParameterPhy> &pars)
{
  const std::string method("std::vector<ParameterPhy*> VantHoff::whatParametersAreWe(const std::vector<ParameterPhy> &)");

  if(pars.size() != 4 && pars.size() != 5)
   antiochError(method,"A Van't Hoff kinetics model needs 4 of 5 parameters, no less, no more...");

  int iE(-1),iD(-1),iB(-1),iT0(-1);
  for(int i = 0; i < (int)pars.size(); i++)
  {
     if(pars[i].isHomogeneous("J/mol"))
     {
        if(iE != -1)iT0 = iE;
        iE = i;
     }
     if(pars[i].unitPar() == "K")
                (iE != -1)?iT0 = i:iE = i;
     if(!pars[i].isUnited())iB = i;
     if(pars[i].unitPar() == "K-1")iD = i;
  }
  if(iT0 == -1 && pars.size() != 4)
        antiochError(method,"It seems you have provided a default parameter but it wasn't found as such...");

  if(iE == -1 || iB == -1 || iD == -1)
    antiochError(method,"At least one the following parameters have not been found: Ea, D or beta.");

  int iA = (((int)pars.size() - 1) * (int)pars.size())/2 - iE - iB - iD;
  if(iT0 != -1)iA -= iT0;
  
  std::vector<const ParameterPhy*> out;
  out.push_back(&pars[iE]);
  out.push_back(&pars[iD]);
  out.push_back(&pars[iB]);
  (iT0 == -1)?out.push_back(&PhyCon::TrefKin):out.push_back(&pars[iT0]);
  out.push_back(&pars[iA]);

  return out;
}

VantHoff &VantHoff::operator=(const VantHoff &rhs)
{
  if(this == &rhs)return *this;
  setPreExp(rhs.getPreExp());
  setBeta(rhs.getBeta());
  setEa(rhs.getEa());
  setD(rhs.getD());
//
  setT0(rhs.getT0());
//
  setT(rhs.getT());

  return *this;
}

template<class M>
M VantHoff::kEquation(const M& temp, const M& pre, const M& beta, const M& ea, const M& d, const M& Tref) const
{
  using std::pow;
  using Antioch::pow;
  using std::exp;
  using Antioch::exp;
  return pre * pow(temp/Tref,beta) * exp(-ea/temp + d * temp);
}
}
