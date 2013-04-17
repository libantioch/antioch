#include "antioch/DuplicateProcess.hpp"

namespace Antioch{
DuplicateProcess::DuplicateProcess(const DuplicateProcess &rhs)
{
  setKineticsProcess("Duplicate process");
  init(rhs.getKineticsModel(),rhs.getAllParameters());
  setTemperature(rhs.getTemperature());
}

DuplicateProcess::DuplicateProcess(const std::string &kinMod, int nPro):
  Process(kinMod)
{
  setKineticsProcess("Duplicate process");
  const std::string method("DuplicateProcess::DuplicateProcess(const std::string &, int ):");
  for(int i = 1; i < nPro; i++)
  {
     rates.push_back(kineticsModeling::modelChoice(kinMod));
  }
  if((int)rates.size() != nPro)antiochError(method,"Problem in constructor!!!");
}

DuplicateProcess::DuplicateProcess(const std::string &kinMod, const std::vector<std::vector<ParameterPhy> > &pars)
{
  setKineticsProcess("Duplicate process");
  init(kinMod,pars);
}

DuplicateProcess::~DuplicateProcess()
{
  for(unsigned int i = 0; i < rates.size(); i++)
  {
     delete rates[i];
     rates[i] = NULL;
  }
  rates.clear();
}

void DuplicateProcess::init(const std::string &kinMod, const std::vector<std::vector<ParameterPhy> > &pars)
{
  setKineticsModel(kinMod);
  setAllParameters(pars);
}

const std::vector<std::vector<ParameterPhy> > DuplicateProcess::getAllParameters() const
{
  std::vector<std::vector<ParameterPhy> > out;
  for(unsigned int i = 0; i < rates.size(); i++)
  {
     out.push_back(getParametersFromRate(rates[i]));
  }
  return out;
}

const std::vector<ParameterPhy> DuplicateProcess::getParameters() const
{
  std::vector<ParameterPhy> out(getParametersFromRate(rates[0]));
  for(unsigned int i = 1; i < rates.size(); i++)
  {
    std::vector<ParameterPhy> tmp(getParametersFromRate(rates[i]));
    for(unsigned int j = 0; j < tmp.size(); j++)
    {
      out.push_back(tmp[j]);
    }
  }
  return out;
}

void DuplicateProcess::setAllParameters(const std::vector<std::vector<ParameterPhy> > &pars)
{
  rates.clear();
  for(unsigned int i = 0; i < pars.size(); i++)
  {
     rates.push_back(kineticsModeling::modelChoice(kineticsMod,pars[i]));
  }
}

void DuplicateProcess::setTemperature(ParameterPhy *T, int nT)
{
  const std::string method("void DuplicateProcess::setTemperature(ParameterPhy *)");
  if(rates[0] == NULL)
        antiochError(method,"Get a rate constant model before you set the temperature!!!");

  for(unsigned int i = 0; i < rates.size(); i++)
  {
    rates[i]->setT(T,nT);
  }
}

double DuplicateProcess::getRateConstantT(double T, int i) const
{
  const std::string method("void DuplicateProcess::setTemperature(ParameterPhy *)");
  if(rates[0] == NULL)
        antiochError(method,"Get a rate constant model before you set the temperature!!!");

  double rateC = 0.;
  for(unsigned int j = 0; j < rates.size(); j++)
  {
      rateC += rates[j]->getRateConstantT(T,i);
  }

  return rateC;
}

std::vector<ParameterPhy> DuplicateProcess::allRateConstant() const
{
  const std::string method("std::vector<ParameterPhy> DuplicateProcess::allRateConstant() const");
  if(rates[0] == NULL)
        antiochError(method,"Get a rate constant model before you calculate the temperature!!!");

  std::vector<ParameterPhy> out;

  for(int i = 0; i < rates[0]->getPreExp().nValues(); i++)out.push_back(rateConstant(i));

  return out;
}

ParameterPhy DuplicateProcess::rateConstant(int nk) const
{
  const std::string method("ParameterPhy DuplicateProcess::rateConstant() const");
  if(rates[0] == NULL)
        antiochError(method,"Get a rate constant model before you calculate the temperature!!!");
  ParameterPhy rateC(rates[0]->rateConstant(nk));
  for(unsigned int i = 1; i < rates.size(); i++)
  {
      rateC += rates[i]->rateConstant(nk);
  }

  return rateC;
}

/*!\todo Rethink it
 *
 */
std::vector<std::vector<ParameterPhy> > DuplicateProcess::sortParameters(const std::vector<ParameterPhy> &pars)
{
  const std::string method("std::vector<std::vector<ParameterPhy> > DuplicateProcess::sortParameters(const std::vector<ParameterPhy> &)");

//A are only parameters having s-1
//A test is (storage->at(i).getUnitObject()->getSIPower("s") == -1)
  std::vector<std::vector<ParameterPhy> > sort;
  sort.resize(1);
  for(unsigned int i = 0; i < pars.size(); i++)
  {
     Units test(pars[i].unitPar());
     if(test.getSIPower("s") == -1)
     {
        sort[0].push_back(pars[i]); //we store the As, the size gives the degree of the duplication
     }
  }

  if((int)pars.size()%sort[0].size() != 0)
    antiochError(method,"Give exactly the same number of parameters for each process please, fix on the todo list");

  std::vector<ParameterPhy> tmp;
  bool next(false);
  for(unsigned int j = 0; j < pars.size(); j++)//look for pairing unit
  {
    next = false;
    tmp.clear();
    for(unsigned int i = 0; i < sort.size(); i++)//compare to what we already have
    {
      if(sort[i][0].unitPar() == pars[j].unitPar())//if the same, out
      {
         next=true;
         break;
      }
    }
    if(next)continue;
    tmp.push_back(pars[j]);//found a new
    for(unsigned int k = j + 1; k < pars.size(); k++)//where are the other(s)
    {
       if(pars[k].unitPar() == tmp[0].unitPar())tmp.push_back(pars[k]);
    }
    if(tmp.size() != sort[0].size())antiochError(method,"Damn you, the algorithm is not good enough!!!");
    sort.push_back(tmp); //add the nProcesses new parameters
  }
  std::vector<std::vector<ParameterPhy> > par;
  par.resize(sort[0].size());
  for(unsigned int i = 0; i < sort.size(); i++) //reverse the order sort[i][j] = par[j][i]
  {
     for(unsigned int j = 0; j < sort[i].size(); j++)
     {
       par[j].resize(sort.size());
       par[j][i].replace(sort[i][j]);
     }
  }

  return par;
}
}
