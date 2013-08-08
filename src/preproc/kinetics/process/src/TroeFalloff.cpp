#include "antioch/TroeFalloff.hpp"

namespace Antioch{
double TroeFalloff::d(0.14);

TroeFalloff::TroeFalloff(const TroeFalloff &rhs)
{
  setKineticsProcess("Troe falloff");
  init(rhs.getKineticsModel(),rhs.getParametersFromRate(rates[0]),rhs.getParametersFromRate(rates[1]));
  setTemperature(rhs.getTemperature());
  setConcentration(rhs.getConcentration());
  setAlpha(rhs.getAlpha());
  setT1(rhs.getT1());
  setT2(rhs.getT2());
  setT3(rhs.getT3());
}

TroeFalloff::TroeFalloff(const std::string &kinMod, 
                         const std::vector<ParameterPhy> &pars0, 
                         const std::vector<ParameterPhy> &parsinf, 
                         const ParameterPhy &al,
                         const ParameterPhy &To,
                         const ParameterPhy &Tt,
                         const ParameterPhy &Tth,
                         ParameterPhy *m):
  FalloffProcess(kinMod,pars0,parsinf,m),
  alpha(al),
  T1(To),
  T2(Tt),
  T3(Tth),
  alphaCal(al),
  T1Cal(To),
  T2Cal(Tt),
  T3Cal(Tth)
{setKineticsProcess("Troe falloff");}

void TroeFalloff::setTemperature(ParameterPhy *Tptr)
{
  rates[0]->setT(Tptr);
  rates[1]->setT(Tptr);
  T = Tptr;
}

const std::vector<ParameterPhy> TroeFalloff::getParameters() const
{
  std::vector<ParameterPhy> out(getParametersFromRate(rates[0]));
  std::vector<ParameterPhy> tmp(getParametersFromRate(rates[1]));
  for(unsigned int j = 0; j < tmp.size(); j++)
  {
    out.push_back(tmp[j]);
  }

  out.push_back(alpha);
  out.push_back(T1);
  out.push_back(T2);
  out.push_back(T3);

  return out;
}

void TroeFalloff::setAlpha(const ParameterPhy &al)
{
  alpha = al;
  alpha.setNamePar("alpha");
  alphaCal = alpha;
}

void TroeFalloff::setT1(const ParameterPhy &t)
{
  T1 = t;
  T1.setNamePar("T*");
  T1Cal = T1;
}
void TroeFalloff::setT2(const ParameterPhy &t)
{
  T2 = t;
  T2.setNamePar("T**");
  T2Cal = T2;
}

void TroeFalloff::setT3(const ParameterPhy &t)
{
  T3 = t;
  T3.setNamePar("T***");
  T3Cal = T3;
}

void TroeFalloff::setFalloffParameters(const std::vector<ParameterPhy> &pars)
{
  const std::string method("void TroeFalloff::setFalloffParameters(const std::vector<ParameterPhy> &)");
  if(pars.size() != 3 && pars.size() != 4)
     antiochError(method,"Please provide 3 (alpha, T*** and T*) or 4 (T**) parameters");

  int iA(-1),iT1(-1),iT2(-1),iT3(-1);
  for(unsigned int i = 0; i < pars.size(); i++)
  {
     if(pars[i].namePar() == "alpha")iA=i;
     if(pars[i].namePar() == "T*")iT1=i;
     if(pars[i].namePar() == "T**")iT2=i;
     if(pars[i].namePar() == "T***")iT3=i;
  }
  if(iA  == -1)antiochError(method,"alpha parameter not found");
  if(iT1 == -1)antiochError(method,"T* parameter not found");
  if(iT3 == -1)antiochError(method,"T*** parameter not found");
  setAlpha(pars[iA]);
  setT3(pars[iT3]);
  setT1(pars[iT1]);
  if(iT2 != -1)setT2(pars[iT2]);
}

void TroeFalloff::resetF(std::vector<double> &pars)
{
  const std::string method("");
  if(pars.size() < nFPars() ||
     pars.size() > nFTotalPars())
                antiochError(method,"Not the correct number of parameters");

  resetAlpha(pars[0]);
  resetT3(pars[1]);
  resetT1(pars[2]);
  if(pars.size() == nFTotalPars())resetT2(pars[3]);
}

const ParameterPhy TroeFalloff::getDefaultAlphaParameter() const
{
  ParameterPhy def("Alpha default parameter",getDefaultAlphaMin(),0.,CORE_UNCERTAINTY_TYPE_NONE,getDefaultAlphaUnit());
  def.addValue(getDefaultAlphaMax());

  return def;
}

const ParameterPhy TroeFalloff::getDefaultT1Parameter() const
{
  ParameterPhy def("T* default parameter",getDefaultT1Min(),0.,CORE_UNCERTAINTY_TYPE_NONE,getDefaultT1Unit());
  def.addValue(getDefaultT1Max());

  return def;
}

const ParameterPhy TroeFalloff::getDefaultT2Parameter() const
{
  ParameterPhy def("T** default parameter",getDefaultT2Min(),0.,CORE_UNCERTAINTY_TYPE_NONE,getDefaultT2Unit());
  def.addValue(getDefaultT2Max());

  return def;
}

const ParameterPhy TroeFalloff::getDefaultT3Parameter() const
{
  ParameterPhy def("T*** default parameter",getDefaultT3Min(),0.,CORE_UNCERTAINTY_TYPE_NONE,getDefaultT3Unit());
  def.addValue(getDefaultT3Max());

  return def;
}

const std::string TroeFalloff::getDefaultFUnit(int i) const 
{
  const std::string method("const std::string TroeFalloff::getDefaultFUnit(int ) const ");
  switch(i)
  {
     case 0:
     {
       return getDefaultAlphaUnit();
       break;
     }
     case 1:
     {
       return getDefaultT3Unit();
       break;
     }
     case 2:
     {
       return getDefaultT1Unit();
       break;
     }
     case 3:
     {
       return getDefaultT2Unit();
       break;
     }
  }
  antiochError(method,"Asking for a non existing parameter");

  return std::string();
}

double TroeFalloff::getDefaultFMin(int i) const
{
  const std::string method("double TroeFalloff::getDefaultFMin(int ) const");
  switch(i)
  {
     case 0:
     {
       return getDefaultAlphaMin();
       break;
     }
     case 1:
     {
       return getDefaultT3Min();
       break;
     }
     case 2:
     {
       return getDefaultT1Min();
       break;
     }
     case 3:
     {
       return getDefaultT2Min();
       break;
     }
  }
  antiochError(method,"Asking for a non existing parameter");

  return 0.;
}

double TroeFalloff::getDefaultFMax(int i) const
{
  const std::string method("double TroeFalloff::getDefaultFMax(int ) const");
  switch(i)
  {
     case 0:
     {
       return getDefaultAlphaMax();
       break;
     }
     case 1:
     {
       return getDefaultT3Max();
       break;
     }
     case 2:
     {
       return getDefaultT1Max();
       break;
     }
     case 3:
     {
       return getDefaultT2Max();
       break;
     }
  }
  antiochError(method,"Asking for a non existing parameter");

  return 0.;
}

const ParameterPhy TroeFalloff::getDefaultFPriorParameter(int i) const
{
  const std::string method("const ParameterPhy TroeFalloff::getDefaultFPriorParameter(int ) const");
  switch(i)
  {
     case 0:
     {
       return getDefaultAlphaParameter();
       break;
     }
     case 1:
     {
       return getDefaultT3Parameter();
       break;
     }
     case 2:
     {
       return getDefaultT1Parameter();
       break;
     }
     case 3:
     {
       return getDefaultT2Parameter();
       break;
     }
  }
  antiochError(method,"Asking for a non existing parameter");

  return ParameterPhy();
}

double TroeFalloff::getDefaultPriorMinValue(int i) const
{
  const std::string method("double TroeFalloff::getDefaultPriorMinValue(int ) const");
  if(rates[0] == NULL)antiochError(method,"No default configuration on non existing rate constant.");

  if(i < (int)nKinPars())return rates[0]->getDefaultMin(i);
  return getDefaultFMin(i-nKinPars());
}

double TroeFalloff::getDefaultPriorMaxValue(int i) const
{
  const std::string method("double TroeFalloff::getDefaultPriorMaxValue(int ) const");
  if(rates[0] == NULL)antiochError(method,"No default configuration on non existing rate constant.");

  if(i < (int)nKinPars())return rates[0]->getDefaultMax(i);
  return getDefaultFMax(i-nKinPars());
}

const std::string TroeFalloff::getDefaultPriorUnit(int i) const
{
  const std::string method("const std::string TroeFalloff::getDefaultPriorUnit(int ) const");
  if(rates[0] == NULL)antiochError(method,"No default configuration on non existing rate constant.");

  if(i < (int)nKinPars())return rates[0]->getDefaultUnit(i);
  return getDefaultFUnit(i-nKinPars());
}

const ParameterPhy TroeFalloff::getDefaultPriorParameter(int i) const
{
  const std::string method("const ParameterPhy TroeFalloff::getDefaultPriorParameter(int ) const");
  if(rates[0] == NULL)antiochError(method,"No default configuration on non existing rate constant.");

  if(i < (int)nKinPars())return rates[0]->getDefaultParameter(i);
  return getDefaultFPriorParameter(i-nKinPars());
}

//logF
ParameterPhy TroeFalloff::logFT(double t, int i)
{
  return logFEquation<ParameterPhy>(ParameterPhy(logFcentT(t,i)),PrT(t,i),
                                    ParameterPhy(cT(t,i)),ParameterPhy(nT(t,i)));
}

ParameterPhy TroeFalloff::logFM(double m, int i)
{
  return logFEquation<ParameterPhy>(logFcent(i),PrM(m,i),c(i),n(i));
}

double TroeFalloff::logFTM(double t, double m, int i) const
{
  return logFEquation<double>(logFcentT(t,i),PrTM(t,m,i),cT(t,i),nT(t,i));
}

//Fcent
ParameterPhy TroeFalloff::logFcent(int nF)
{
  const std::string method("ParameterPhy TroeFalloff::Fcent()");
  if(alpha.empty())antiochError(method,alpha.namePar() + " is not defined!");
  if(T1.empty())antiochError(method,T1.namePar() + " is not defined!");
  if(T3.empty())antiochError(method,T3.namePar() + " is not defined!");
  if(T == NULL)antiochError(method,"The temperature is not defined!");

  resetPar(alphaCal,alpha.value(nF));
  resetPar(T1Cal,T1.value(nF));
  if(!T2.empty())resetPar(T2Cal,T2.value(nF));
  resetPar(T3Cal,T3.value(nF));

  ParameterPhy out = logFcentEquation<ParameterPhy>(*T, T3Cal, T1Cal, alphaCal,T2Cal);

  if(alpha.isUncDef())out.correcUncAddBiais(2.0 * exp(-(*T)/T3-(*T)/T1),alpha);

  return out;
}

double TroeFalloff::logFcentT(double t, int i) const
{
  const std::string method("double TroeFalloff::FcentT(double , int )");
  if(alpha.empty())antiochError(method,alpha.namePar() + " is not defined!");
  if(T1.empty())antiochError(method,T1.namePar() + " is not defined!");
  if(T3.empty())antiochError(method,T3.namePar() + " is not defined!");

  double t2 = 0.;
  if(!T2.empty())t2 = T2.value(i);
  
  return logFcentEquation<double>(t, T3.value(i), T1.value(i), alpha.value(i),t2);
}

template <class P>
P TroeFalloff::logFEquation(const P& logFc, const P& pr, const P& C, const P& N) const
{
  using std::log10;
  using Antioch::log10;
  using std::pow;
  using Antioch::pow;
  return ( logFc /
             (
                P(1.) + pow( (log10(pr) + C) / (N - (log10(pr) + C ) * d) , 2)
             )
         );
}

template <class P>
P TroeFalloff::CEquation(const P& Fc) const
{
   using std::log10;
   using Antioch::log10;
   return (P(-0.4) - P(0.67) * log10(Fc));
}

template <class P>
P TroeFalloff::NEquation(const P& logFc) const
{
   return (P(0.75) - P(1.27) * logFc);
}

template <class P>
P TroeFalloff::logFcentEquation(const P& temp, const P& t3, const P& t1, const P& al, const P& t2) const
{
   using std::log10;
   using Antioch::log10;
   using std::exp;
   using Antioch::exp;
   P out = exp(-temp/t3) * (P(1.0) - al) +  exp(-temp/t1) * al;
   if(!T2.empty()) out += exp(-t2/temp);
   return log10(out);
}
}
