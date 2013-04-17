#include "antioch/KineticsModel.hpp"

namespace Antioch{
void KineticsModel::checkA(const std::string &object) const
{
  checkParam(A,"A",object);
}

void KineticsModel::checkT(const std::string &object) const
{
  const std::string method("void KineticsModel::checkT(const std::string &)");
  if(T == NULL)antiochError(method,"T is undefined in kinetics model " + object + "!!!");
}

void KineticsModel::checkParam(const ParameterPhy &test, const std::string &namepar, const std::string &object) const
{
  if(test.empty())antiochError(object,namepar + " has no values");
}

void KineticsModel::setT(ParameterPhy *Tused, int nT)
{
  const std::string method("void Reaction::setT(ParameterPhy *, int (= 0))");
  if(nT != 0)antiochError(method,"You need a multiple Hercour-Essen for that... (nT != 0)");

  T = Tused;
}

ParameterPhy * KineticsModel::getT(int nT) const
{
  const std::string method("ParameterPhy * Reaction::getT(int  (= 0)) const");
  if(nT != 0)antiochError(method,"You need a multiple Hercour-Essen for that... (nT != 0)");

  return T;
}

const std::vector<ParameterPhy> KineticsModel::allRateConstant()
{
  std::vector<ParameterPhy> out;

  for(int i = 0; i < A.nValues(); i++)
  {
    out.push_back(rateConstant(i));
  }

  return out;
}

std::vector<const ParameterPhy*> KineticsModel::findSameUnitInTwoParameters(const std::vector<ParameterPhy> &pars,const std::string &unit)
{
  const std::string method("std::vector<ParameterPhy*> KineticsModel::findUnitInTwoParameters(const std::vector<ParameterPhy> &,const std::string &)");
  if(pars.size() != 2)antiochError(method,"Only two parameters for this method please.");
  std::vector<const ParameterPhy*> out;
  if(pars[0].unitPar() == unit)
  {
      out.push_back(&pars[0]);
      out.push_back(&pars[1]);
  }else if(pars[1].unitPar() == unit)
  {
      out.push_back(&pars[1]);
      out.push_back(&pars[0]);
  }

  return out;
}

std::vector<const ParameterPhy*> KineticsModel::findHomoUnitInTwoParameters(const std::vector<ParameterPhy> &pars,const std::string &unit)
{
  const std::string method("std::vector<ParameterPhy*> KineticsModel::findHomoUnitInTwoParameters(const std::vector<ParameterPhy> &,const std::string &)");
  if(pars.size() != 2)antiochError(method,"Only two parameters for this method please.");
  std::vector<const ParameterPhy*> out;
  if(pars[0].isHomogeneous(unit))
  {
      out.push_back(&pars[0]);
      out.push_back(&pars[1]);
  }else if(pars[1].isHomogeneous(unit))
  {
      out.push_back(&pars[1]);
      out.push_back(&pars[0]);
  }

  return out;
}

const ParameterPhy KineticsModel::getDefaultPreExp() const
{
  ParameterPhy def("Pre-exponential parameter",getDefaultPreExpMin(),0.,CORE_UNCERTAINTY_TYPE_NONE,getDefaultPreExpUnit());
  def.addValue(getDefaultPreExpMax());

  return def;
}

void KineticsModel::resetPar(ParameterPhy &par, double r)
{
  par.clearValues();
  par.addValue(r);
}

void KineticsModel::showAll(std::ostream &out) const
{
  out << "*****************" << std::endl;
  showKinModel(out);
  showT(out);
  out << "*****************" << std::endl;
}

void KineticsModel::showT(std::ostream &out) const
{
  if(T ==NULL)
  {
     out << "No temperature set" << std::endl;
  }else
  {
     out << "Temperature parameter at adress: " << T << std::endl;
     T->showAll(out);
  }
}
}
