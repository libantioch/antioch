#include "antioch/Arrhenius.hpp"
#include "antioch/phymath.hpp"
#include "antioch/physical_constants.hpp"

namespace Antioch{

Arrhenius::Arrhenius(const std::vector<ParameterPhy> &pars)
{
  init(pars);
}

void Arrhenius::init(const std::vector<ParameterPhy> &pars)
{
  std::vector<const ParameterPhy*> sortedPar(whatParametersAreWe(pars));

  setEa(*sortedPar[0]);  
  setPreExp(*sortedPar[1]);
}

void Arrhenius::reset(const std::vector<double> &pars)
{
  const std::string method("void Arrhenius::reset(const std::vector<double> &)");

//number of parameters
  if(pars.size() != nPars())antiochError(method,"Not the right number of parameters");
  resetPreExp(pars[0]);
  resetEa(pars[1]);
}

void Arrhenius::resetEa(double ae)
{
  resetPar(Ea,ae);
  setExpPa();
}

const std::vector<ParameterPhy> Arrhenius::getParameters() const
{
  std::vector<ParameterPhy> out;
  out.push_back(A);
  out.push_back(Ea);

  return out;
}

const ParameterPhy Arrhenius::getDefaultEa() const
{
  ParameterPhy def("Default Ea parameter",getDefaultEaMin(),0.,CORE_UNCERTAINTY_TYPE_NONE,getDefaultEaUnit());
  def.addValue(getDefaultEaMax());

  return def;
}

const std::string Arrhenius::getDefaultUnit(int i) const
{
  const std::string method("const std::string Arrhenius::getDefaultUnit(int ) const");
  switch(i)
  {
     case 0:
     {
       return getDefaultPreExpUnit();
       break;
     }
     case 1:
     {
       return getDefaultEaUnit();
       break;
     }
  }
  antiochError(method,"Asking for a non existing parameter");

  return std::string();
}

double Arrhenius::getDefaultMin(int i) const
{
  const std::string method("double Arrhenius::getDefaultMin(int ) const");
  switch(i)
  {
     case 0:
     {
       return getDefaultPreExpMin();
       break;
     }
     case 1:
     {
       return getDefaultEaMin();
       break;
     }
  }
  antiochError(method,"Asking for a non existing parameter");

  return 0.;
}

double Arrhenius::getDefaultMax(int i) const
{
  const std::string method("double Arrhenius::getDefaultMax(int ) const");
  switch(i)
  {
     case 0:
     {
       return getDefaultPreExpMax();
       break;
     }
     case 1:
     {
       return getDefaultEaMax();
       break;
     }
  }
  antiochError(method,"Asking for a non existing parameter");

  return 0.;
}

const ParameterPhy Arrhenius::getDefaultParameter(int i) const
{
  const std::string method("const ParameterPhy Arrhenius::getDefaultParameter(int )");
  switch(i)
  {
     case 0:
     {
       return getDefaultPreExp();
       break;
     }
     case 1:
     {
       return getDefaultEa();
       break;
     }
  }
  antiochError(method,"Asking for a non existing parameter");

  return ParameterPhy();
}

std::vector<const ParameterPhy*> Arrhenius::whatParametersAreWe(const std::vector<ParameterPhy> &pars)
{
  const std::string method("std::vector<ParameterPhy*> Arrhenius::whatParametersAreWe(const std::vector<ParameterPhy> &)");
  if(pars.size() != 2)
      antiochError(method,"Arrhenius-type of kinetics rate constant use 2 parameters: A and Ea. Please provide exactly two parameters for this kinetics model.");
  std::vector<const ParameterPhy*> out;

  out = findSameUnitInTwoParameters(pars,"K");

  if(out.size() == 0)
  {
    out = findHomoUnitInTwoParameters(pars,"J/mol");
  }
  if(out.size() != 2)
      antiochError(method,"No activation energy has been found (either K united or J/mol homogeneous), what is this for an Arrhenius?");

  return out;
}

void Arrhenius::setEa(const ParameterPhy &ae)
{
  Ea.replace(ae);
  setExpPa();
}

void Arrhenius::checkParams() const
{
  checkA("Arrhenius::checkParams() const");
  checkParam(Ea,"Ea","Arrhenius::checkParams() const");
}

const ParameterPhy Arrhenius::rateConstant(int nk)
{
  checkT("Arrhenius");
  checkParams();
  resetPar(ACal,A.value(nk));
  resetPar(expPaCal,expPa.value(nk));
  return kEquation<ParameterPhy>(*T,ACal,expPaCal);
}

void Arrhenius::setExpPa()
{
  const std::string method("void Arrhenius::setExpPa()");
  expPa = Ea;
  if(expPa.isHomogeneous("J/mol"))expPa /= PhyCon::R;
  if(!expPa.isHomogeneous("K"))antiochError(method,"Problem in the activation energy. Its unit is " 
                                                 + Ea.unitPar() 
                                                 + ". It should be either an energy per mole (something homogeneous to J/mol) or K.");
  expPa.changeToSomeUnit("K");
  expPaCal = expPa;
}

double Arrhenius::getRateConstantT(double Tem, int i) const
{
  return kEquation<double>(Tem,A.value(i),expPa.value(i));
}

Arrhenius &Arrhenius::operator=(const Arrhenius &rhs)
{
  if(this == &rhs)return *this;
  setPreExp(rhs.getPreExp());
  setEa(rhs.getEa());
//
  setT(rhs.getT());

  return *this;
}

template<class M>
M Arrhenius::kEquation(const M& temp, const M& pre, const M& ea) const
{
  using std::exp;
  using Antioch::exp;
  return  pre * exp(-ea/temp);
}
}
