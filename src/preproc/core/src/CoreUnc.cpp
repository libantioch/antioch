//-----------------------------------------------------------------------bl-
//

//Antioch
#include "antioch/Error.hpp"
#include "antioch/CoreUnc.hpp"

//C++

namespace Antioch{
CoreUnc::CoreUnc(double dval, std::string typeStr):
   type(typeStr)
{
  if(type == CORE_UNCERTAINTY_TYPE_ABSOLUTE)
     {
        addDvalue(dval*dval);
     }
  else
     {
        addDvalue(dval);
     }
}

double CoreUnc::getDvalue(int i) const
{
  if(type == CORE_UNCERTAINTY_TYPE_NONE)
  {
    return 0.;
  }else if(type == CORE_UNCERTAINTY_TYPE_ABSOLUTE || type == CORE_UNCERTAINTY_TYPE_RELATIVE)
  {
    if(i >= (int)dvalues.size() && dvalues.size() != 1)
      antiochError("double CoreUnc::getDValue(int) const",
                   "You ask for a value that is not stored for type " + type + ", i > sizeDVal() - 1.");
    return (i < (int)dvalues.size())?dvalues[i]:dvalues[0];
  }else
  {
    antiochError("double CoreUnc::getDValue(int) const",
                   "What is this type \"" + type + "\" ?");
    return -1.;
  }
}

double CoreUnc::getSqrtvalue(int i) const
{
  return std::sqrt(dvalues[i]);
}

void CoreUnc::setDvalue(double dval, int i)
{
  if(dvalues.size() == 0 && i != 0)
     antiochError("void CoreUnc::setDvalue(double,int)",
                   "Trying to assign non first value to empty std::vector.");

  if((int)dvalues.size() == i)
  {
    dvalues.push_back(std::fabs(dval));
  }else
  {
    dvalues[i] = std::fabs(dval);
  }
}

void CoreUnc::addDvalues(const std::vector<double> &dval)
{
  for(unsigned int i = 0; i < dval.size(); i++)
  {
    addDvalue(dval[i]);
  }
}

const std::vector<double> CoreUnc::getSqrtvalues() const
{
  std::vector<double> out;
  for(int i = 0; i < sizeDVal(); i++)
  {
      out.push_back(std::sqrt(dvalues[i]));
  }

  return out;
}

void CoreUnc::setToNullUnc()
{
  for(unsigned int i = 0; i < dvalues.size(); i++)
  {
     dvalues[i] = 0.;
  }
}

void CoreUnc::setTypeToUnknown()
{
  type = CORE_UNCERTAINTY_TYPE_UNKNOWN;
  dvalues.clear();
}

void CoreUnc::setTypeToNone()
{
  type = CORE_UNCERTAINTY_TYPE_NONE;
  dvalues.clear();
  dvalues.push_back(0.);
}

bool const CoreUnc::operator==(const CoreUnc &rhs) const
{
  bool isEqual = ((rhs.sizeDVal() == (int)dvalues.size()) && (type == rhs.getType()));
  if(isEqual)
  {
    for(int i = 0; i < rhs.sizeDVal(); i++)
    {
      isEqual = isEqual && (dvalues[i] == rhs.getDvalue(i));
    }
  }
  return isEqual;
}

CoreUnc & CoreUnc::operator =(const CoreUnc &rhs)
{
  if(this == &rhs)return *this;
  type = rhs.getType();
  dvalues = rhs.getDvalues();
  return *this;
}

CoreUnc & CoreUnc::operator+=(const CoreUnc &rhs)  
{//dealing with errors: types not consistent
  const std::string method("CoreUnc & CoreUnc::operator+=(const CoreUnc &)  ");

//if unknown, it is contagious
  if(rhs.getType() == CORE_UNCERTAINTY_TYPE_UNKNOWN)
  {
    setTypeToUnknown();
    return *this;
  }
//if none or we are unknown, nothing to do
  if(rhs.getType() == CORE_UNCERTAINTY_TYPE_NONE ||
     type == CORE_UNCERTAINTY_TYPE_UNKNOWN)return *this;

//if we are none, just copy his
  if(type == CORE_UNCERTAINTY_TYPE_NONE)
  {
     *this = rhs;
     return *this;
  }

//if relative, it's not possible
  if(type == CORE_UNCERTAINTY_TYPE_RELATIVE ||
     rhs.getType() == CORE_UNCERTAINTY_TYPE_RELATIVE)
  {
    std::string errMess = "There is a type problem.\n";
    errMess += "You want me to add " + type + " and " + rhs.getType();
    errMess += ".\nI cannot add CoreUnc in such conditions.\n";
    antiochError(method,errMess);
  }
// here the type is ok, now we need to test the number of values

// only values possible are
//   - either equal or
//   - only one value for one of them

  int min = (sizeDVal() < rhs.sizeDVal())?sizeDVal():rhs.sizeDVal();
  int max = (sizeDVal() > rhs.sizeDVal())?sizeDVal():rhs.sizeDVal();

  if((min != 1) && (min != max))
  {
    std::string errMess =  "There are not the same values of uncertainties. ";
    errMess += "It is not permitted to add uncertainties this way ";
    errMess += "in CoreUnc objects.\n";
    antiochError(method,errMess);
//now compute only if rhs actually have something to change
  }

  int i,z(0);
  int *k(&i);
  if(sizeDVal() < max)dvalues.resize(max,dvalues[0]);
  if(rhs.sizeDVal() < max)k = &z;
  for(i = 0; i < max; i++)
  {
    dvalues[i] += rhs.getDvalue(*k);
  }

  return *this;
}

CoreUnc & CoreUnc::operator-=(const CoreUnc &rhs)
{
   return (*this += rhs);
}

CoreUnc & CoreUnc::operator*=(double r)
{
   if(type != CORE_UNCERTAINTY_TYPE_ABSOLUTE)return *this;
   for(int i = 0; i < sizeDVal(); i++)
   {
     dvalues[i] *= r*r; 
   }
   return *this;
}

CoreUnc & CoreUnc::operator/=(double r)
{
   return (*this *= (1./r));
}

}
