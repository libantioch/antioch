#include "antioch/BerthelotHercourtEssen.hpp"
#include "antioch/phymath.hpp"

namespace Antioch{
BerthelotHercourtEssen::BerthelotHercourtEssen(const std::vector<ParameterPhy> &pars)
{
  init(pars);
}

void BerthelotHercourtEssen::init(const std::vector<ParameterPhy> &pars)
{
  std::vector<const ParameterPhy*> out = whatParametersAreWe(pars);

  setD(*out[0]);
  setBeta(*out[1]);
  setT0(*out[2]);
  setPreExp(*out[3]);
}

void BerthelotHercourtEssen::reset(const std::vector<double> &pars)
{
  const std::string method("void BerthelotHercourtEssen::reset(const std::vector<double> &)");

//number of parameters
  if(pars.size() != nPars())antiochError(method,"Not the right number of parameters");
  resetPreExp(pars[0]);
  resetBeta(pars[1]);
  resetD(pars[2]);
}

void BerthelotHercourtEssen::checkParams() const
{
  checkA("BerthelotHercourt-Essen::checkParams() const");
  checkParam(beta,"beta","BerthelotHercourtEssen::checkParams() const");
  checkParam(D,"D","BerthelotHercourtEssen::checkParams() const");
}

const ParameterPhy BerthelotHercourtEssen::rateConstant(int nk)
{
  checkT("Berthelot Hercourt-Essen");
  checkParams();
  ParameterPhy out = kEquation<ParameterPhy>(*T,A.value(nk),beta.value(nk),D.value(nk),T0.value(nk));
  out.setNamePar("Berthelot Hercourt-Essen rate constant");
  if(
     (T->isUncDef()) && // T is uncertain
     (beta.value(nk) != 0.)    && // beta is not null
     (D.value(nk) != 0.)          // D is not null
     )
        out.correcUncAddBiais(2.*D.value(nk)*beta.value(nk)/(*T)*pow(out,2),*T);
  return out;
}

double BerthelotHercourtEssen::getRateConstantT(double Tem, int i) const
{
  return kEquation<double>(Tem, A.value(i), beta.value(i), D.value(i), T0.value());
}

const std::vector<ParameterPhy> BerthelotHercourtEssen::getParameters() const
{
  std::vector<ParameterPhy> out;
  out.push_back(A);
  out.push_back(beta);
  out.push_back(D);
  out.push_back(T0);

  return out;
}

const std::string BerthelotHercourtEssen::getDefaultUnit(int i) const
{
  const std::string method("const std::string BerthelotHercourtEssen::getDefaultUnit(int ) const");
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
       return getDefaultDUnit();
       break;
     }
  }
  antiochError(method,"Asking for a non existing parameter");

  return std::string();
}

double BerthelotHercourtEssen::getDefaultMin(int i) const
{
  const std::string method("double BerthelotHercourtEssen::getDefaultMin(int ) const");
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
       return getDefaultDMin();
       break;
     }
  }
  antiochError(method,"Asking for a non existing parameter");

  return 0.;
}

double BerthelotHercourtEssen::getDefaultMax(int i) const
{
  const std::string method("double BerthelotHercourtEssen::getDefaultMax(int ) const");
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
       return getDefaultDMax();
       break;
     }
  }
  antiochError(method,"Asking for a non existing parameter");

  return 0.;
}

const ParameterPhy BerthelotHercourtEssen::getDefaultParameter(int i) const
{
  const std::string method("const ParameterPhy BerthelotHercourtEssen::getDefaultParameter(int )");
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
       return getDefaultD();
       break;
     }
  }
  antiochError(method,"Asking for a non existing parameter");

  return ParameterPhy();
}

std::vector<const ParameterPhy*> BerthelotHercourtEssen::whatParametersAreWe(const std::vector<ParameterPhy> &pars)
{
  const std::string method("std::vector<ParameterPhy*> BerthelotHercourtEssen::whatParametersAreWe(const std::vector<ParameterPhy> &)");

  if(pars.size() < 3 || pars.size() > 4)
        antiochError(method,"Number of parameters does not match Berthelot Hercourt-Essen definition (3 or 4).");

  int iB(-1),iD(-1),iT0(-1);
  for(int i = 0; i < (int)pars.size(); i++)
  {
     if(!pars[i].isUnited())iB = i;
     if(pars[i].unitPar() == "K")iT0 = i;
     if(pars[i].unitPar() == "K-1")iD = i;
  }
  if(iD == -1 || iB == -1)
        antiochError(method,"D (K-1) or beta (no unit) was no found in the parameters.");
  if(iT0 == -1 && pars.size() != 3)
        antiochError(method,"It seems you have provided a default parameter but it wasn't found as such...");

  int iA = ( ((int)pars.size() - 1) * (int)pars.size() )/2 - iB + iD;
  if(iT0 != -1) iA -= iT0;

  std::vector<const ParameterPhy*> out;
  out.push_back(&pars[iD]);
  out.push_back(&pars[iB]);
  (iT0 == -1)?out.push_back(&PhyCon::TrefKin):out.push_back(&pars[iT0]);
  out.push_back(&pars[iA]);

  return out;
}

BerthelotHercourtEssen &BerthelotHercourtEssen::operator=(const BerthelotHercourtEssen &rhs)
{
  if(this == &rhs)return *this;
  setPreExp(rhs.getPreExp());
  setBeta(rhs.getBeta());
  setD(rhs.getD());
//
  setT0(rhs.getT0());
//
  setT(rhs.getT());

  return *this;
}

template<class M>
M BerthelotHercourtEssen::kEquation(const M& temp , const M& pre, const M& beta, const M& d, const M& Tref) const
{
  using std::pow;
  using Antioch::pow;
  using std::exp;
  using Antioch::exp;
  return pre * pow(temp/Tref,beta) * exp(d * temp);
}
}
