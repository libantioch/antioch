#include "antioch/FalloffProcess.hpp"

namespace Antioch{
FalloffProcess::FalloffProcess(const std::string &kinMod, 
                               const std::vector<ParameterPhy> &pars0, 
                               const std::vector<ParameterPhy> &parsinf, 
                               ParameterPhy *m):
 Process(kinMod,pars0)
{
  setConcentration(m);
  rates.push_back(NULL);
  setParameters(parsinf,rates[1]);
}

FalloffProcess::~FalloffProcess()
{
  for(unsigned int i = 0; i < rates.size(); i++)
  {
     delete rates[i];
     rates[i] = NULL;
  }
  if(M != NULL)delete M;
  M = NULL;
}

void FalloffProcess::resetParameters(const std::vector<double> &pars)
{
  const std::string method("void FalloffProcess::resetParameters(const std::vector<double> &)");
  if(rates[0] == NULL)antiochError(method,"No rate constant at lower pressures defined");
  if(rates[1] == NULL)antiochError(method,"No rate constant at higher pressures defined");

  

  if(((int)pars.size()-nFPars())%2 != 0)antiochError(method,"You don't have the correct number of parameters!");

  std::vector<double> pars0,parsinf,parsF;
  for(unsigned int i = 0; i < (pars.size() - (unsigned int)nFPars())/2; i++)
  {
    pars0.push_back(pars[i]);
    parsinf.push_back(pars[i + (pars.size() + (unsigned int)nFPars())/2]);
  }
  for(int i = ((int)pars.size() - nFPars()); i < (int)pars.size(); i++)
  {
    parsF.push_back(pars[i]);
  }

  rates[0]->reset(pars0);
  rates[1]->reset(parsinf);
  resetF(parsF);
}

const std::vector<ParameterPhy> FalloffProcess::getParameters() const
{
  std::vector<ParameterPhy> out(getParametersFromRate(rates[0]));
  std::vector<ParameterPhy> tmp(getParametersFromRate(rates[1]));
  for(unsigned int j = 0; j < tmp.size(); j++)
  {
    out.push_back(tmp[j]);
  }

  return out;
}

double FalloffProcess::PrTM(double t, double m, int i) const
{
  const std::string method("double FalloffProcess::Pr(double , double , int )");
  if(rates[0] == NULL)
        antiochError(method,"Low-pressure rate constant k0 undefined!!");
  if(rates[1] == NULL)
        antiochError(method,"High-pressure rate constant kinf undefined!!");

  return PrEquation<double>(rates[0]->getRateConstantT(t,i),rates[1]->getRateConstantT(t,i),m);
}

ParameterPhy FalloffProcess::PrT(double t, int i) const
{
  const std::string method("double FalloffProcess::Pr(double , int )");
  if(rates[0] == NULL)
          antiochError(method,"Low-pressure rate constant k0 undefined!!");
  if(rates[1] == NULL)
          antiochError(method,"High-pressure rate constant kinf undefined!!");
  if(M == NULL)
          antiochError(method,"Concentration of mixture undefined!!");

  return PrEquation<ParameterPhy>(ParameterPhy("k0",  rates[0]->getRateConstantT(t,i),0.,CORE_UNCERTAINTY_TYPE_NONE,rates[0]->getPreExp().unitPar()),
                                  ParameterPhy("kinf",rates[1]->getRateConstantT(t,i),0.,CORE_UNCERTAINTY_TYPE_NONE,rates[1]->getPreExp().unitPar()),
                                  *M);
}

ParameterPhy FalloffProcess::PrM(double m, int i) const
{
  const std::string method("double FalloffProcess::Pr(double , int )");
  if(rates[0] == NULL)
        antiochError(method,"Low-pressure rate constant k0 undefined!!");
  if(rates[1] == NULL)
        antiochError(method,"High-pressure rate constant kinf undefined!!");

  return PrEquation<ParameterPhy>(rates[0]->rateConstant(i),rates[1]->rateConstant(i),
                                  ParameterPhy("M",m,0.,CORE_UNCERTAINTY_TYPE_NONE,defaultMUnit()));
}

void FalloffProcess::setTemperature(ParameterPhy *t)
{
  rates[0]->setT(t);
  rates[1]->setT(t);
}

void FalloffProcess::init(const std::string &kinMod, const std::vector<ParameterPhy> &pars0, const std::vector<ParameterPhy> &parsinf, ParameterPhy *m)
{
  Process::init(kinMod,pars0); //so we know what we're talking about
  setParameters(parsinf,rates[1]);
  M = m;
}

ParameterPhy FalloffProcess::rateConstantT(double t, int i)
{
  const std::string method("double FalloffProcess::getRateConstantT(double , int) const");
  if(rates[0] == NULL)
        antiochError(method,"Low-pressure rate constant k0 undefined!!");
  if(rates[1] == NULL)
        antiochError(method,"High-pressure rate constant kinf undefined!!");
  if(M == NULL)
        antiochError(method,"Mixture concentration not defined!!");
  
  return kEquation<ParameterPhy>(ParameterPhy("kinf",rates[1]->getRateConstantT(t,i),0.,CORE_UNCERTAINTY_TYPE_NONE,rates[0]->getPreExp().unitPar()),
                                 PrT(t,i),
                                 exp(logFT(t,i)));
}

ParameterPhy FalloffProcess::rateConstantM(double m, int i)
{
  const std::string method("double FalloffProcess::getRateConstantM(double , int) const");
  if(rates[0] == NULL)
        antiochError(method,"Low-pressure rate constant k0 undefined!!");
  if(rates[1] == NULL)
        antiochError(method,"High-pressure rate constant kinf undefined!!");
  
  return kEquation<ParameterPhy>(rates[1]->rateConstant(i),PrM(m,i),exp(logFM(m,i)));
}


double FalloffProcess::getRateConstantTM(double t, double m, int i) const
{
  const std::string method("double FalloffProcess::getRateConstantTM(double , double, int) const");
  if(rates[0] == NULL)
        antiochError(method,"Low-pressure rate constant k0 undefined!!");
  if(rates[1] == NULL)
        antiochError(method,"High-pressure rate constant kinf undefined!!");
  
  return kEquation<double>(rates[1]->getRateConstantT(t,i),PrTM(t,m,i),std::exp(logFTM(t,m,i)));
}

template<class N>
N FalloffProcess::PrEquation(const N &lowpres, const N &highpres, const N& conc) const
{
  return conc * lowpres / highpres;
}

template<class N>
N FalloffProcess::kEquation(const N &highpres, const N &pr, const N& f) const
{
  return highpres * pr / (1.0 * pr) * f;
}
}
