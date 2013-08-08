#include "antioch/HercourtEssen.hpp"
#include "antioch/phymath.hpp"

namespace Antioch{

HercourtEssen::HercourtEssen(const std::vector<ParameterPhy> &pars)
{
  init(pars);
}

void HercourtEssen::init(const std::vector<ParameterPhy> &pars)
{
  std::vector<const ParameterPhy*> sortedPar(whatParametersAreWe(pars));

  setBeta(*sortedPar[0]);
  setT0(*sortedPar[1]);
  setPreExp(*sortedPar[2]);
}

void HercourtEssen::reset(const std::vector<double> &pars)
{
  const std::string method("void HercourtEssen::reset(const std::vector<double> &)");

//number of parameters
  if(pars.size() != nPars())antiochError(method,"Not the right number of parameters");
  resetPreExp(pars[0]);
  resetBeta(pars[1]);
}

void HercourtEssen::checkParams() const
{
  checkA("Hercourt-Essen::checkParams() const");
  checkParam(beta,"beta","Hercourt-Essen::checkParams() const");
}

const ParameterPhy HercourtEssen::rateConstant(int nk)
{
  checkT("Hercourt-Essen");
  checkParams();
  resetPar(ACal,A.value(nk));
  resetPar(betaCal,beta.value(nk));
  return kEquation<ParameterPhy>(*T, ACal, betaCal, T0);
}

double HercourtEssen::getRateConstantT(double Tem, int i) const
{
  return kEquation<double>(Tem,A.value(i), beta.value(i), T0.value());
}

const std::vector<ParameterPhy> HercourtEssen::getParameters() const
{
  std::vector<ParameterPhy> out;
  out.push_back(A);
  out.push_back(beta);
  out.push_back(T0);

  return out;
}

const ParameterPhy HercourtEssen::getDefaultBeta() const
{
  ParameterPhy def("Default beta configuration",getDefaultBetaMin(),0.,CORE_UNCERTAINTY_TYPE_NONE,getDefaultBetaUnit());
  def.addValue(getDefaultBetaMax());
  return def;
}

const std::string HercourtEssen::getDefaultUnit(int i) const
{
  const std::string method("const std::string HercourtEssen::getDefaultUnit(int ) const");
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
  }
  antiochError(method,"Asking for a non existing parameter");

  return std::string();
}

double HercourtEssen::getDefaultMin(int i) const
{
  const std::string method("double getDefaultMin(int ) const");
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
  }
  antiochError(method,"Asking for a non existing parameter");

  return 0.;
}

double HercourtEssen::getDefaultMax(int i) const
{
  const std::string method("double getDefaultMax(int ) const");
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
  }
  antiochError(method,"Asking for a non existing parameter");

  return 0.;
}

const ParameterPhy HercourtEssen::getDefaultParameter(int i) const
{
  const std::string method("const ParameterPhy HercourtEssen::getDefaultParameter(int )");
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
  }
  antiochError(method,"Asking for a non existing parameter");

  return ParameterPhy();
}


std::vector<const ParameterPhy*> HercourtEssen::whatParametersAreWe(const std::vector<ParameterPhy> &pars)
{
  const std::string method("std::vector<ParameterPhy*> HercourtEssen::whatParametersAreWe(const std::vector<ParameterPhy> &)");
  if(pars.size() != 2 && pars.size() != 3)
      antiochError(method,"Hercourt-Essen type of kinetics rate constant use 2 or 3 parameters: A, beta and possibly T0, if not 300 K. Please provide exactly two parameters for this kinetics model.");

  int iT0(-1),iB(-1);
  for(unsigned int i = 0; i < pars.size(); i++)
  {
     if(pars[i].unitPar() == "K")iT0 = (int)i;
     if(!pars[i].isUnited())iB = (int)i;
  }
  if(iB == -1)
      antiochError(method,"No power parameter has been found (no unit), what is this for a Hercourt-Essen?");

  if(iT0 == -1 && pars.size() == 3)
        antiochError(method,"It seems you have provided a default parameter but it wasn't found as such...");

  int iA = (((int)pars.size() - 1) * (int)pars.size())/2 - iB;
  if(iT0 != -1)iA -= iT0;
     
  std::vector<const ParameterPhy*> out;
  out.push_back(&pars[iB]);
  (iT0 == -1)?out.push_back(&PhyCon::TrefKin):out.push_back(&pars[iT0]);
  out.push_back(&pars[iA]);

  return out;
}

void HercourtEssen::showKinModel(std::ostream &out) const
{
  out << "Hercourt-Essen model" << std::endl;
  A.showAll(out);
  beta.showAll(out);
  T0.showAll(out);

}

HercourtEssen &HercourtEssen::operator=(const HercourtEssen &rhs)
{
  if(this == &rhs)return *this;
  setPreExp(rhs.getPreExp());
  setBeta(rhs.getBeta());
//
  setT0(rhs.getT0());
//
  setT(rhs.getT());

  return *this;
}

template<class M>
M HercourtEssen::kEquation(const M& temp, const M& pre, const M& p, const M& tref) const
{
  using std::pow;
  using Antioch::pow;
  return pre * pow(temp/tref,p);
}
}
