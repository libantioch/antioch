#include "antioch/MultHercourtEssen.hpp"
#include "antioch/phymath.hpp"

namespace Antioch{

MultHercourtEssen::MultHercourtEssen(const std::vector<ParameterPhy> &pars)
{
   init(pars);
}

void MultHercourtEssen::init(const std::vector<ParameterPhy> &pars)
{
   std::vector<const ParameterPhy *> spars = whatParametersAreWe(pars);

   setBetae(*spars[0]);
   setBetai(*spars[1]);
   setT0e(*spars[2]);
   setT0i(*spars[3]);
   setPreExp(*spars[4]);
}

void MultHercourtEssen::checkParams() const
{
  const std::string method("void MultHercourtEssen::checkParams() const");
  checkA("MultHercourt-Essen::checkParams() const");
  checkParam(betae,"betae","MultHercourt-Essen::checkParams() const");
  checkParam(betai,"betae","MultHercourt-Essen::checkParams() const");
  if(Ti == NULL)antiochError(method,"The second temperature is undefined.");
}

void MultHercourtEssen::reset(const std::vector<double> &pars)
{
  const std::string method("void MultHercourtEssen::reset(const std::vector<double> &)");

//number of parameters
  if(pars.size() != nPars())antiochError(method,"Not the right number of parameters");
  resetPreExp(pars[0]);
  resetBetae(pars[1]);
  resetBetai(pars[2]);
}


void MultHercourtEssen::setT(ParameterPhy *Tused, int nT)
{
  const std::string method("void MultHercourtEssen::setT(ParameterPhy *, int (=0))");
  if(nT > 1)antiochError(method,"Multiple Hercourt-Essen has (yet) only two temperatures (nT <= 1)!!");

  if(nT == 0)
  {
     setTe(Tused);
  }else
  {
     setTi(Tused);
  }
}

ParameterPhy * MultHercourtEssen::getT(int nT) const
{
  const std::string method("ParameterPhy * MultHercourtEssen::getT(int (= 0)) const");
  if(nT > 1)antiochError(method,"Multiple Hercourt-Essen has (yet) only two temperatures (nT <= 1)!!");

  if(nT == 0)
  {
     return getTe();
  }else
  {
     return getTi();
  }
}

const ParameterPhy MultHercourtEssen::getDefaultBetae() const
{
  ParameterPhy def("Default betae configuration",getDefaultBetaeMin(),0.,CORE_UNCERTAINTY_TYPE_NONE,getDefaultBetaeUnit());
  def.addValue(getDefaultBetaeMax());
  return def;
}

const ParameterPhy MultHercourtEssen::getDefaultBetai() const
{
  ParameterPhy def("Default betai configuration",getDefaultBetaiMin(),0.,CORE_UNCERTAINTY_TYPE_NONE,getDefaultBetaiUnit());
  def.addValue(getDefaultBetaiMax());
  return def;
}

const ParameterPhy MultHercourtEssen::rateConstant(int nk)
{
  checkT("multiple Hercourt-Essen");
  checkParams();
  resetPar(ACal,A.value(nk));
  resetPar(betaeCal,betae.value(nk));
  resetPar(betaiCal,betai.value(nk));
  return kEquation<ParameterPhy>(*T, *Ti, ACal, betaeCal, betaiCal, T0e, T0i);
}

double MultHercourtEssen::getRateConstantT(double Tem, int i) const
{
  antiochError("double MultHercourtEssen::getRateConstantT(double , int ) const",
            "This method should not be used. Use a Hercourt-Essen model if you want only one temperature.");
  return -1.;
}

double MultHercourtEssen::getRateConstantT(double Tem, double Tim, int i) const
{
  return kEquation<double>(Tem, Tim, A.value(i),betae.value(i),betai.value(i), T0e.value(), T0i.value());
}

const std::vector<ParameterPhy> MultHercourtEssen::getParameters() const
{
  std::vector<ParameterPhy> out;
  out.push_back(A);
  out.push_back(betae);
  out.push_back(betai);

  return out;
}

const std::string MultHercourtEssen::getDefaultUnit(int i) const
{
  const std::string method("const std::string MultHercourtEssen::getDefaultUnit(int ) const");
  switch(i)
  {
     case 0:
     {
       return getDefaultPreExpUnit();
       break;
     }
     case 1:
     {
       return getDefaultBetaeUnit();
       break;
     }
     case 2:
     {
       return getDefaultBetaiUnit();
       break;
     }
  }
  antiochError(method,"Asking for a non existing parameter");

  return std::string();
}

double MultHercourtEssen::getDefaultMin(int i) const
{
  const std::string method("double MultHercourtEssen::getDefaultMin(int ) const");
  switch(i)
  {
     case 0:
     {
       return getDefaultPreExpMin();
       break;
     }
     case 1:
     {
       return getDefaultBetaeMin();
       break;
     }
     case 2:
     {
       return getDefaultBetaiMin();
       break;
     }
  }
  antiochError(method,"Asking for a non existing parameter");

  return 0.;
}

double MultHercourtEssen::getDefaultMax(int i) const
{
  const std::string method("double MultHercourtEssen::getDefaultMax(int ) const");
  switch(i)
  {
     case 0:
     {
       return getDefaultPreExpMax();
       break;
     }
     case 1:
     {
       return getDefaultBetaeMax();
       break;
     }
     case 2:
     {
       return getDefaultBetaiMax();
       break;
     }
  }
  antiochError(method,"Asking for a non existing parameter");

  return 0.;
}

const ParameterPhy MultHercourtEssen::getDefaultParameter(int i)
{
  const std::string method("const ParameterPhy MultHercourtEssen::getDefaultParameter(int )");
  switch(i)
  {
     case 0:
     {
       return getDefaultPreExp();
       break;
     }
     case 1:
     {
       return getDefaultBetae();
       break;
     }
     case 2:
     {
       return getDefaultBetai();
       break;
     }
  }
  antiochError(method,"Asking for a non existing parameter");

  return ParameterPhy();
}

std::vector<const ParameterPhy*> MultHercourtEssen::whatParametersAreWe(const std::vector<ParameterPhy> &pars)
{
  const std::string method("std::vector<ParameterPhy*> MultHercourtEssen::whatParametersAreWe(const std::vector<ParameterPhy> &)");
  if(pars.size() < 3 || pars.size() > 5)
     antiochError(method,"You need to provide at least 3 parameters and at most 5 for a multiple Hercourt-Essen kinetics model.");

  int iB0(-1),iB1(-1),iT00(-1),iT01(-1);
  for(int i = 0; i < (int)pars.size(); i++)
  {
    if(!pars[i].isUnited() && iB0 == -1)iB0 = i;
    if(!pars[i].isUnited() && iB0 != -1)iB1 = i;
    if(pars[i].unitPar() == "K" && iT01 == -1)iT00 = i;
    if(pars[i].unitPar() == "K" && iT01 != -1)iT01 = i;
  }
  if(iB0 == -1 || iB1 == -1)
     antiochError(method,"At least one power parameter was not found.");

  if(iT00 == -1 && pars.size() != 3)
        antiochError(method,"It seems you have provided a default parameter but it wasn't found as such...");
  if(iT01 == -1 && iT00 != -1 && pars.size() != 4)
        antiochError(method,"It seems you have provided a default parameter but it wasn't found as such...");

  int iA = ( ((int)pars.size() -1) * (int)pars.size() )/2 - iB0 - iB1; //n(n+1)/2, sum of all the indexes, then we substract the found
  
  if(iT01 != -1)
  { 
     iA -= iT00 + iT01;
  }else if(iT00 != -1)
  {
     iA -= iT00;
  }

  std::vector<const ParameterPhy*> out;
  out.push_back(&pars[iB0]);
  out.push_back(&pars[iB1]);
  (iT00 == -1)?out.push_back(&PhyCon::TrefKin):out.push_back(&pars[iT00]);
  (iT01 == -1)?out.push_back(out.back()):out.push_back(&pars[iT01]);
  out.push_back(&pars[iA]);

  return out;
}

MultHercourtEssen &MultHercourtEssen::operator=(const MultHercourtEssen &rhs)
{
  if(this == &rhs)return *this;
  setPreExp(rhs.getPreExp());
  setBetae(rhs.getBetae());
  setBetai(rhs.getBetai());
//
  setT0e(rhs.getT0e());
  setT0i(rhs.getT0i());
//
  setT(rhs.getT());
  setTi(rhs.getTi());

  return *this;
}

template<class M>
M MultHercourtEssen::kEquation(const M& temp1, const M& temp2, const M& pre, const M& beta1, const M& beta2, const M& Tref1, const M& Tref2) const
{
  using std::pow;
  using Antioch::pow;
  return (pre * pow(temp1/Tref1,beta1) * pow(temp2/Tref2,beta2));
}
}
