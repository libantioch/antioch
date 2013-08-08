#include "antioch/Process.hpp"
namespace Antioch{

Process::Process(const Process &rhs):
  M(NULL)
{
  rates.push_back(NULL);
  init(rhs.getKineticsModel(),rhs.getParameters());
  setTemperature(rhs.getTemperature());
}

Process::Process(const std::string &kinMod, const std::vector<ParameterPhy> &pars):
  M(NULL)
{
  rates.push_back(NULL);
  init(kinMod,pars);
}

Process::~Process()
{
  if(rates[0] != NULL)delete rates[0];
}

void Process::init(const std::string &kinMod, const std::vector<ParameterPhy> &pars)
{
  setKineticsModel(kinMod);
  setParameters(pars,rates[0]);
}

const std::vector<ParameterPhy> Process::getParametersFromRate(KineticsModel * rateC) const
{
  const std::string method("const std::vector<ParameterPhy> Process::getParametersFromRate(KineticsModel * )");
  if(rateC == NULL)antiochError(method,"What parameters can you possibly want from an unexisting rate constant?");
  return rateC->getParameters();
}

std::vector<ParameterPhy> Process::allRateConstant() const
{
  const std::string method("std::vector<ParameterPhy> allRateConstant() const");
  if(rates[0] == NULL)
        antiochError(method,"Get a rate constant model before you calculate its values!!!");

  return rates[0]->allRateConstant();
}

double Process::getRateConstantTM(double t, double m, int i) const
{
  const std::string method("ParameterPhy Process::getRateConstantTM(double, double , int (=0)) const");

  antiochError(method,"What are you doing here??");

  return 0.;
}

ParameterPhy Process::rateConstantM(double m, int i) const
{
  const std::string method("ParameterPhy Process::getRateConstantM(double , int (=0)) const");

  antiochError(method,"What are you doing here??");

  return ParameterPhy(0.);
}

ParameterPhy Process::rateConstantT(double t, int i) const
{
  const std::string method("ParameterPhy Process::getRateConstantT(double , int (=0)) const");

  antiochError(method,"What are you doing here??");

  return ParameterPhy(0.);
}

ParameterPhy Process::rateConstant(int nk) const
{
  const std::string method("ParameterPhy Process::rateConstant(int (= 0)) const;");
  if(rates[0] == NULL)
        antiochError(method,"Get a rate constant model before you calculate its values!!!");

  return rates[0]->rateConstant(nk);
}

double Process::getRateConstantT(double T, int i) const
{
  const std::string method("double Process::getRateConstantT(double , int ) const;");
  if(rates[0] == NULL)
        antiochError(method,"Get a rate constant model before a value!!!");

  return rates[0]->getRateConstantT(T,i);
}

void Process::setTemperature(ParameterPhy *t, int nT)
{
  const std::string method("void Process::setTemperature(ParameterPhy *, int ( = 0))");
  if(rates[0] == NULL)
        antiochError(method,"Get a rate constant model before you set the temperature!!!");

  rates[0]->setT(t,nT);
}

ParameterPhy *Process::getTemperature(int nT) const
{
  const std::string method("ParameterPhy *Process::getTemperature()");
  if(rates[0] == NULL)
        antiochError(method,"Get a rate constant model before you want to have the temperature!!!");

  return rates[0]->getT(nT);
}

bool Process::isTDefined()
{
  const std::string method("ParameterPhy Process::isTDefined()");
  if(rates[0] == NULL)
        antiochError(method,"Get a rate constant model before you want to test the temperature!!!");

  return rates[0]->isTDefined();
}

void Process::setParameters(const std::vector<ParameterPhy> &pars, KineticsModel *&rateC)
{
  if(rateC == NULL)
  {
     rateC = kineticsModeling::modelChoice(kineticsMod,pars);
  }else
  {
    rateC->init(pars);
  }
}

void Process::resetPar(ParameterPhy &par,double d)
{
  par.clearParameter();
  par.addValue(d);
}

void Process::resetParameters(const std::vector<double> &pars)
{
  const std::string method("void Process::resetParameters(const std::vector<double> &)");

  int nparPerproc(pars.size()/getNProcesses()),k(0);
  for(unsigned int iproc = 0; iproc < getNProcesses(); iproc++)
  {
    std::vector<double> parsTmp;
    for(int i = 0; i < nparPerproc; i++)
    {
       parsTmp.push_back(pars[k]);
       k++;
    }
    if(rates[iproc] == NULL)antiochError(method,"No rate constant defined");
    rates[iproc]->reset(parsTmp);
  }
}

unsigned int Process::nKinPars(int iProc) const
{
  const std::string method("unsigned int Process::nKinPars(int (= 0)) const");
  if(iProc >= (int)getNProcesses())antiochError(method,"Asking for non existing process (iProc > getNProcesses())");
  if(rates[iProc] == NULL)antiochError(method,"Rate constant undefined");

  return rates[iProc]->nPars();
}

unsigned int Process::nKinOptionPars(int iProc) const
{
  const std::string method("unsigned int Process::nKinOptionPars(int (= 0)) const");
  if(iProc >= (int)getNProcesses())antiochError(method,"Asking for non existing process (iProc > getNProcesses())");
  if(rates[iProc] == NULL)antiochError(method,"Rate constant undefined");

  return rates[iProc]->nOptionPars();
}

unsigned int Process::nTotalKinPars(int iProc) const
{
  const std::string method("unsigned int Process::nTotalKinPars(int (= 0)) const");
  if(iProc >= (int)getNProcesses())antiochError(method,"Asking for non existing process (iProc > getNProcesses())");
  if(rates[iProc] == NULL)antiochError(method,"Rate constant undefined");

  return rates[iProc]->nTotalPars();
}

double Process::getDefaultPriorMinValue(int i) const
{
  const std::string method("double Process::defaultPriorMinValue(int ) const"); 
  if(rates[0] == NULL)antiochError(method,"No default configuration on non existing rate constant.");

  return rates[0]->getDefaultMin(i);
}

double Process::getDefaultPriorMaxValue(int i) const
{
  const std::string method("double Process::defaultPriorMaxValue(int ) const"); 
  if(rates[0] == NULL)antiochError(method,"No default configuration on non existing rate constant.");

  return rates[0]->getDefaultMax(i);
}

const std::string Process::getDefaultPriorUnit(int i) const
{
  const std::string method("const std::string Process::defaultPriorUnit(int ) const"); 
  if(rates[0] == NULL)antiochError(method,"No default configuration on non existing rate constant.");

  return rates[0]->getDefaultUnit(i);
}

const ParameterPhy Process::getDefaultPriorParameter(int i) const
{
  const std::string method("const ParameterPhy Process::defaultPriorParameter(int ) const"); 
  if(rates[0] == NULL)antiochError(method,"No default configuration on non existing rate constant.");

  return rates[0]->getDefaultParameter(i);
}

void Process::showAll(std::ostream &out) const
{
  out << "Chemical process: " << kineticsProc << std::endl;
  out << "Kinetics model: "   << kineticsMod  << std::endl;
  if(M == NULL)
  {
    out << "No pressure set"<< std::endl;
  }else
  {
     out << "Pressure parameter at adress: "  << M << std::endl;
     M->showAll(out);
  }

  for(unsigned int i = 0; i < rates.size(); i++)
  {
     if(rates[i] == NULL)
     {
        out << "No rate constant" << std::endl;
        continue;
     }
     rates[i]->showAll(out);
  }
}
}
