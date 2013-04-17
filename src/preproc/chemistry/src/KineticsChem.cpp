#include "antioch/KineticsChem.hpp"

namespace Antioch{

KineticsChem::~KineticsChem()
{
  if(reactionRate != NULL)delete reactionRate;
}

KineticsChem & KineticsChem::operator=(const KineticsChem &rhs)
{
  if(this == &rhs)return *this;
  initialize(rhs.getSimpleChemPart(),rhs.getKineticsProcess(),rhs.getReaction());

  return *this;
}

void KineticsChem::initialize(const SimpleChem &SC, const std::string &kinProc, const Process * kinObj)
{
  const std::string method("void KineticsChem::initialize(const SimpleChem &, const std::string &, const Process * )");
  setNonReactivePart(SC);
  setKineticsProcess(kinProc);
  if(kinObj != NULL)
  {
    kinModel = kinObj->getKineticsModel();
    if(kinProcess != kinObj->getKineticsProcess())antiochError(method,"You're trying to set a " + kinProcess + "with a " + kinObj->getKineticsProcess());
    setKinetics(kinObj);
  }
}

void KineticsChem::setKineticsOnly(const KineticsChem &KC)
{
  setKineticsProcess(KC.getKineticsProcess());
  setMolecules(KC.getMolecules());
  if(KC.getReaction() != NULL)
  {
    kinModel = KC.getKineticsModel();
    setKinetics(KC.getReaction());
  }
}

void KineticsChem::initialize(const SimpleChem &SC, const std::string &kinProc, const std::string &kinMod,
                              const std::vector<ParameterPhy> &pars,const std::vector<ParameterPhy> &parsinf)
{
  setNonReactivePart(SC);
  setKineticsProcess(kinProc);
  setKineticsModel(kinMod);
  setKinetics(pars,parsinf);
}

void KineticsChem::setKinetics(const Process *kinObj)
{
  const std::string method("void KineticsChem::setKinetics(const Process *)");
  kinProcess = kinObj->getKineticsProcess();
  if(kinProcess.empty())
        antiochError(method,"Define a chemical process before setting the kinetics model");
  setReactionRate(kinModel,kinObj);
}

void KineticsChem::setKinetics(const std::vector<ParameterPhy> &pars, const std::vector<ParameterPhy> &parsinf)
{
  const std::string 
        method("void KineticsChem::setKinetics(const std::string &, const std::vector<ParameterPhy> &, const std::vector<ParameterPhy> &)");
  if(kinProcess.empty())
        antiochError(method,"Define a chemical process before setting the kinetics model");
  if(kinModel.empty())kinModel = "Kooij";

  setReactionRate(kinModel);
  if(pars.empty())return;
  if(parsinf.empty())
  {
    reactionRate->init(kinModel,pars);
    return;
  }
  if(!isFalloff())antiochError(method,"This is not a falloff you're requesting, why two kinetics sets of parameters?");
  reactionRate->init(kinModel,pars,parsinf);
}

void KineticsChem::setReactionRate(const std::string &kinMod, const Process * kinPtr)
{
  const std::string method("void KineticsChem::setReactionRate(const std::string & (=\"Kooij\"))");
  if(kinProcess.empty())antiochError(method,"No chemical process set.");
  if(reactionRate != NULL)delete reactionRate;
  reactionRate = NULL;
  reactionRate = BUCKchemicalProcess::getChemicalProcess(kinProcess,kinPtr);
}

const std::string KineticsChem::getKineticsModel() const
{
  if(reactionRate == NULL)
  {
     return kinModel;
  }else
  {
    return reactionRate->getKineticsModel();
  }
  return std::string();
}

const std::vector<ParameterPhy> KineticsChem::getKineticsParameters() const
{
  const std::string method("const std::vector<ParameterPhy> KineticsChem::getKineticsParameters() const;");
  if(reactionRate == NULL)
        antiochError(method,"Cannot provide kinetics parameters from a non-existing reaction rate");

  return reactionRate->getParameters();
}

int KineticsChem::nKinPars() const
{
  const std::string method("int KineticsChem::nKinPars() const");
  if(reactionRate == NULL)
        antiochError(method,"What number of kinetics parmaeters can you want from an unexisting reaction rate??");
  return (int)reactionRate->nKinPars();
}

std::vector<ParameterPhy> KineticsChem::allRateConstant() const
{
  const std::string method("std::vector<ParameterPhy> KineticsChem::allRateConstant() const");
  if(reactionRate == NULL)
        antiochError(method,"What rate constant can you want from an unexisting reaction rate??");
  return reactionRate->allRateConstant();
}

ParameterPhy KineticsChem::rateConstant(int nk) const
{
  const std::string method("ParameterPhy KineticsChem::rateConstant(int ) const");
  if(reactionRate == NULL)
        antiochError(method,"What rates constant can you want from an unexisting reaction rate??");
  return reactionRate->rateConstant(nk);
}

double KineticsChem::getRateConstantT(double t, int i) const
{
  const std::string method("double KineticsChem::getRateConstantT(double ) const");
  if(reactionRate == NULL)
        antiochError(method,"What rate constant can you want from an unexisting reaction rate??");
  return reactionRate->getRateConstantT(t,i);
}

ParameterPhy KineticsChem::rateConstantT(double t, int i)
{
  const std::string method("ParameterPhy KineticsChem::rateConstantT(double , int )");
  if(!isFalloff())antiochError(method,"You need a falloff process for this. You have a " + kinProcess + ".");
  if(reactionRate == NULL)
        antiochError(method,"What rate constant can you want from an unexisting reaction rate??");

  return reactionRate->rateConstantT(t,i);
}

ParameterPhy KineticsChem::rateConstantM(double m, int i)
{
  const std::string method("ParameterPhy KineticsChem::rateConstantM(double , int )");
  if(!isFalloff())antiochError(method,"You need a falloff process for this. You have a " + kinProcess + ".");
  if(reactionRate == NULL)
        antiochError(method,"What rate constant can you want from an unexisting reaction rate??");

  return reactionRate->rateConstantM(m,i);
}

double KineticsChem::getRateConstantTM(double t, double m, int i) const
{
  const std::string method("double KineticsChem::getRateConstantTM(double, double , int ) const");
  if(!isFalloff())antiochError(method,"You need a falloff process for this. You have a " + kinProcess + ".");
  if(reactionRate == NULL)
        antiochError(method,"What rate constant can you want from an unexisting reaction rate??");

  return (reactionRate->getRateConstantTM(t,m,i));
}

bool KineticsChem::isFalloff() const
{
  return (kinProcess.find("falloff") != std::string::npos);
}

void KineticsChem::storeParameter (const ParameterPhy &par)
{
  if(storage == NULL)storage = new std::vector<ParameterPhy>;
  storage->push_back(par);
}

void KineticsChem::storeParameters(const std::vector<ParameterPhy> &pars)
{
  if(storage == NULL)storage = new std::vector<ParameterPhy>;
  for(unsigned int i = 0; i < pars.size(); i++)
  {
      storage->push_back(pars[i]);
  }
}

void KineticsChem::setKineticsTemperature(ParameterPhy *temp)
{
  const std::string method("void KineticsChem::setKineticsTemperature(ParameterPhy *)");
  if(reactionRate == NULL)antiochError(method,"Cannot set the temperature on an unexisting rate constant");

  reactionRate->setTemperature(temp);
}

void KineticsChem::setKineticsConcentration(ParameterPhy *conc)
{
  const std::string method("void KineticsChem::setKineticsConcentration(ParameterPhy *);");
  if(reactionRate == NULL)antiochError(method,"Cannot set the concentration on an unexisting rate constant");

  reactionRate->setConcentration(conc);
}

void KineticsChem::setReactionFromStorage()
{
  const std::string method("void KineticsChem::setReactionFromStorage()");
  if(kinProcess.empty())
     antiochError(method,"No chemical process defined, cannot set a reaction rate.");

  if(storage == NULL)
     antiochError(method,"Storage not defined, cannot set a reaction rate.");


  if(kinModel.empty())
  {
     antiochWarning(method,"No kinetics model stored, Kooij used.");
     kinModel = "Kooij";
  }
  if(reactionRate != NULL)delete reactionRate;
  
  reactionRate = BUCKchemicalProcess::getChemicalProcess(kinProcess);

  if(storage->empty())
  {
     antiochWarning(method,"Chemical process exists but no kinetics parameters!");
     reactionRate->setKineticsModel(kinModel);
     return;
  }

  if(!isFalloff())
  {
    reactionRate->init(kinModel,*storage);
  }else
  {
    std::vector<ParameterPhy> k0,kinf,Fpars;
    parseFalloffParameters(k0,kinf,Fpars);
    reactionRate->init(kinModel,k0,kinf);
    reactionRate->setFalloffParameters(Fpars);
  }

  delete storage;
  storage = NULL;
}

void KineticsChem::parseFalloffParameters(std::vector<ParameterPhy> &k0, std::vector<ParameterPhy> &kinf)
{
  std::vector<ParameterPhy> emp;
  parseFalloffParameters(k0,kinf,emp);
}

void KineticsChem::showKinetics(std::ostream &out) const
{
  const std::string method("void KineticsChem::showKinetics(const std::ostream & (= cout)) const");
  if(reactionRate == NULL)
  {
    out << "No reaction rate" << std::endl;
  }else
  {
    reactionRate->showAll(out);
  }
}

void KineticsChem::parseFalloffParameters(std::vector<ParameterPhy> &k0, std::vector<ParameterPhy> &kinf, std::vector<ParameterPhy> &Fpars)
{
  const std::string method("void KineticsChem::parseFalloffParameters(std::vector<ParameterPhy> &, std::vector<ParameterPhy> &, std::vector<ParameterPhy> &)");

//Troe, parameters are named alpha, T*, T** and T***, imposed by TroeFalloff class
  if(reactionRate->nFPars() != 0)
  {
    unsigned int i(0),j(0),nFfound(0);
    while(i < storage->size() - 1)
    {
      for(i = j; i < storage->size(); i++)
      {
        if(storage->at(i).namePar() == "alpha" ||
           storage->at(i).namePar() == "T*"    ||
           storage->at(i).namePar() == "T**"   ||
           storage->at(i).namePar() == "T***")
        {
          Fpars.push_back(storage->at(i));
          nFfound++;
          storage->erase(storage->begin() + i);
          j = i;
          break;
        }
      }//end of storage loop
    }//end of while loop
    if(nFfound < reactionRate->nFPars() || 
       nFfound > reactionRate->nFTotalPars())antiochError(method,"Haven't found an appropriate number of F parameters.");
  }//end of F parameters


///kinetics part
//k0 and kinf tell apart with A, A in kinf has unit of Ak0 * [M], thus
//Akinf/Ak0 homogeneous to m3 or m3/mol ([M]-1)
//
//Also, A are only parameter having s-1
//A test is (storage->at(i).getUnitObject()->getSIPower("s") == -1)
  int nA(0),iA(-1),iiA(-1);
  for(unsigned int i = 0; i < storage->size(); i++)
  {
     Units test(storage->at(i).unitPar());
     if(test.getSIPower("s") == -1)
     {
        if(nA == 2)antiochError(method,"Too many pre-exponential parameters!");
        nA++;
        (nA == 1)?iA = i:iiA = i;
     }
     if(nA == 2)break;
  }

  bool k0first(false);
  Units tmp(storage->at(iA).unitPar() + "/(" + storage->at(iiA).unitPar() + ")");
  if(tmp.getSIPower("m") < 0)// Ak0 / Akinf
  {
     k0.push_back(storage->at(iA));
     kinf.push_back(storage->at(iiA));
     k0first = true;
  }else if(tmp.getSIPower("m") > 0)// Akinf / Ak0
  {
     k0.push_back(storage->at(iiA));
     kinf.push_back(storage->at(iA));
     k0first = false;
  }else
  {
    antiochError(method,"What are those units for pre-exponential parameters?\n" + 
                      storage->at(iA).unitPar() + "\n" +
                      storage->at(iiA).unitPar());
  }
//the rest is done as following:
//find pairing parameters (same units)
//the found order Ak0 Akinf is applied to
//the other parameters
  unsigned int nFound(2);
  for(unsigned i = 0; i < storage->size(); i++)
  {
     if((int)i == iA || (int)i == iiA)continue;
     for(unsigned int j = i; j < storage->size(); j++)
     {
       if(storage->at(j).isHomogeneous(storage->at(i).unitPar()))//found pair
       {
          nFound++;
          nFound++;
          if(k0first)
          {
            k0.push_back(storage->at(i));
            kinf.push_back(storage->at(j));
          }else
          {
            k0.push_back(storage->at(j));
            kinf.push_back(storage->at(i));
          }
       }
     }
  if(nFound == storage->size())break;
  }

  if(nFound != storage->size())antiochError(method,"It does not work");
}
}
