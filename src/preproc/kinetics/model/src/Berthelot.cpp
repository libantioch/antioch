#include "antioch/Berthelot.hpp"

namespace Antioch{
Berthelot::Berthelot(const std::vector<ParameterPhy> &pars)
{
  init(pars);
}

void Berthelot::init(const std::vector<ParameterPhy> &pars)
{
  std::vector<const ParameterPhy*> out = whatParametersAreWe(pars);

  setD(*out[0]);
  setPreExp(*out[1]);
}

void Berthelot::checkParams() const
{
  checkA("Berthelot::checkParams() const");
  checkParam(D,"D","Berthelot::checkParams() const");
}

void Berthelot::reset(const std::vector<double> &pars)
{
  const std::string method("void Berthelot::reset(const std::vector<double> &)");

//number of parameters
  if(pars.size() != nPars())antiochError(method,"Not the right number of parameters");
  resetPreExp(pars[0]);
  resetD(pars[1]);
}

const ParameterPhy Berthelot::rateConstant(int nk) 
{
  checkT("Berthelot");
  checkParams();
  resetPar(ACal,A.value(nk));
  resetPar(DCal,D.value(nk));
  return kEquation<ParameterPhy>(*T,ACal,DCal);
}

double Berthelot::getRateConstantT(double Tem, int i) const
{
  return kEquation<double>(Tem,A.value(i),D.value(i));
}

const std::vector<ParameterPhy> Berthelot::getParameters() const
{
  std::vector<ParameterPhy> out;
  out.push_back(A);
  out.push_back(D);

  return out;
}

const ParameterPhy Berthelot::getDefaultD() const
{
  ParameterPhy def("Default D parameter",getDefaultDMin(),0.,CORE_UNCERTAINTY_TYPE_NONE,getDefaultDUnit());
  def.addValue(getDefaultDMax());

  return def;
}

const std::string Berthelot::getDefaultUnit(int i) const
{
  const std::string method("const std::string Berthelot::getDefaultUnit(int ) const");
  switch(i)
  {
     case 0:
     {
       return getDefaultPreExpUnit();
       break;
     }
     case 1:
     {
       return getDefaultDUnit();
       break;
     }
  }
  antiochError(method,"Asking for a non existing parameter");

  return std::string();
}

double Berthelot::getDefaultMin(int i) const
{
  const std::string method("double Berthelot::getDefaultMin(int ) const");
  switch(i)
  {
     case 0:
     {
       return getDefaultPreExpMin();
       break;
     }
     case 1:
     {
       return getDefaultDMin();
       break;
     }
  }
  antiochError(method,"Asking for a non existing parameter");

  return 0.;
}

double Berthelot::getDefaultMax(int i) const
{
  const std::string method("double Berthelot::getDefaultMax(int ) const");
  switch(i)
  {
     case 0:
     {
       return getDefaultPreExpMax();
       break;
     }
     case 1:
     {
       return getDefaultDMax();
       break;
     }
  }
  antiochError(method,"Asking for a non existing parameter");

  return 0.;
}

const ParameterPhy Berthelot::getDefaultParameter(int i) const
{
  const std::string method("const ParameterPhy Berthelot::getDefaultParameter(int )");
  switch(i)
  {
     case 0:
     {
       return getDefaultPreExp();
       break;
     }
     case 1:
     {
       return getDefaultD();
       break;
     }
  }
  antiochError(method,"Asking for a non existing parameter");

  return ParameterPhy();
}

std::vector<const ParameterPhy*> Berthelot::whatParametersAreWe(const std::vector<ParameterPhy> &pars)
{
  const std::string method("std::vector<ParameterPhy*> whatParametersAreWe(const std::vector<ParameterPhy> &)");
  if(pars.size() != 2)
        antiochError(method,"A Berthelot model needs 2 and only 2 parameters. Check your entries.");

  std::vector<const ParameterPhy*> out = findSameUnitInTwoParameters(pars,"K-1");

  if(out.size() != 2)
      antiochError(method,"No exponential parameter has been found (K-1 united), what is this for a Berthelot?");

  return out;
}

Berthelot &Berthelot::operator=(const Berthelot &rhs)
{
  if(this == &rhs)return *this;
  setPreExp(rhs.getPreExp());
  setD(rhs.getD());
//
  setT(rhs.getT());

  return *this;
}

template<class M>
M Berthelot::kEquation(const M& temp, const M& pre, const M& d) const
{
  using std::exp;
  using Antioch::exp;
  return pre * exp(temp * d);
}
}
