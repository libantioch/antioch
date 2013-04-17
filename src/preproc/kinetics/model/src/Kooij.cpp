#include "antioch/Kooij.hpp"
#include "antioch/phymath.hpp"

namespace Antioch{

Kooij::Kooij(const std::vector<ParameterPhy> &pars)
{
  init(pars);
}

void Kooij::init(const std::vector<ParameterPhy> &pars)
{
  std::vector<const ParameterPhy*> spars = whatParametersAreWe(pars);

  setEa(*spars[0]);
  setBeta(*spars[1]);
  setT0(*spars[2]);
  setPreExp(*spars[3]);
}

void Kooij::reset(const std::vector<double> &pars)
{
  const std::string method("void Kooij::reset(const std::vector<double> &)");

//number of parameters
  if(pars.size() != nPars())antiochError(method,"Not the right number of parameters");
  resetPreExp(pars[0]);
  resetBeta(pars[1]);
  resetEa(pars[2]);
}

void Kooij::checkParams() const
{
  checkA("Kooij::checkParams");
  checkParam(beta,"beta","Kooij::checkParams() const");
  checkParam(Ea,"Ea","Kooij::checkParams() const");
}

const ParameterPhy Kooij::rateConstant(int nk)
{
  checkT("Kooij");
  checkParams();

  resetPar(ACal,A.value(nk));
  resetPar(betaCal,beta.value(nk));
  resetPar(expPaCal,expPa.value(nk));

  ParameterPhy out = kEquation<ParameterPhy>(*T,ACal,betaCal,expPaCal,T0);
  out.setNamePar("Kooij rate constant");

  if(
     (T->isUncDef())                              && //T is uncertain
     (beta.value(nk) != 0.) && // beta is not null
     (Ea.value(nk)   != 0.)    // Ea is not null
    )
      out.correcUncAddBiais(2.*pow(out,2)*betaCal*expPaCal/(pow(*T,3)),*T);

  return out;
}

double Kooij::getRateConstantT(double Tem, int i) const
{
  return kEquation<double>(Tem,A.value(i),beta.value(i),expPa.value(i),T0.value());
}

const std::vector<ParameterPhy> Kooij::getParameters() const
{
  std::vector<ParameterPhy> out;
  out.push_back(A);
  out.push_back(beta);
  out.push_back(Ea);
  out.push_back(T0);

  return out;
}

const std::string Kooij::getDefaultUnit(int i) const
{
  const std::string method("const std::string Kooij::getDefaultUnit(int ) const");
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
  }
  antiochError(method,"Asking for a non existing parameter");

  return std::string();
}

double Kooij::getDefaultMin(int i) const
{
  const std::string method("double Kooij::getDefaultMin(int ) const");
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
  }
  antiochError(method,"Asking for a non existing parameter");

  return 0.;
}

double Kooij::getDefaultMax(int i) const
{
  const std::string method("double Kooij::getDefaultMax(int ) const");
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
  }
  antiochError(method,"Asking for a non existing parameter");

  return 0.;
}

const ParameterPhy Kooij::getDefaultParameter(int i) const
{
  const std::string method("const Kooij::ParameterPhy getDefaultParameter(int )");
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
  }
  antiochError(method,"Asking for a non existing parameter");

  return ParameterPhy();
}

std::vector<const ParameterPhy*> Kooij::whatParametersAreWe(const std::vector<ParameterPhy> &pars)
{
  const std::string method("std::vector<ParameterPhy*> Kooij::whatParametersAreWe(const std::vector<ParameterPhy> &)");
  if(pars.size() != 3 && pars.size() != 4)
   antiochError(method,"A Kooij kinetics model require 3 or 4 parameters.");

  int iE(-1),iT0(-1),iB(-1);
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
  }
  if(iB == -1 || iE == -1)
        antiochError(method,"Ea or beta (or both) was (were) not found.");
  if(iT0 == -1 && pars.size() != 3)
        antiochError(method,"It seems you have provided a default parameter but it wasn't found as such...");
  int iA = ( ((int)pars.size() - 1) * (int)pars.size() )/2 - iE - iB;
  if(iT0 != -1)iA -= iT0;

  std::vector<const ParameterPhy*> out;
  out.push_back(&pars[iE]);
  out.push_back(&pars[iB]);
  (iT0 == -1)?out.push_back(&PhyCon::TrefKin):out.push_back(&pars[iT0]);
  out.push_back(&pars[iA]);

  return out;
}

Kooij &Kooij::operator=(const Kooij &rhs)
{
  if(this == &rhs)return *this;
  setPreExp(rhs.getPreExp());
  setBeta(rhs.getBeta());
  setEa(rhs.getEa());
//
  setT0(rhs.getT0());
//
  setT(rhs.getT());

  return *this;
}

void Kooij::showKinModel(std::ostream &out) const
{
  out << "Kooij model" << std::endl;
  A.showAll(out);
  beta.showAll(out);
  Ea.showAll(out);
  T0.showAll(out);
}

template<class M>
M Kooij::kEquation(const M& temp, const M& pre, const M& beta, const M& ea, const M& tref) const
{
  using std::pow;
  using Antioch::pow;
  using std::exp;
  using Antioch::exp;
  return pre * pow(temp/tref,beta) * exp(-ea/temp);
}
}
