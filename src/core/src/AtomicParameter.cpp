//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

//Antioch
#include "antioch/Error.hpp"
#include "antioch/AtomicParameter.hpp"

//C++


namespace Antioch{

const std::string AtomicParameter::computationError(std::string ope,std::string op) const
{
  std::string str =
    "I do not allow " + ope + "ing AtomicParameter when there is "
 +  "not the same number of values, and I have not found a "
 +  "constant parameter (only one value). You are trying to " + ope + " "
 +  "parameters with different values.size(), this is not "
 +  "permissible in AtomicParameter. If you want it to be, "
 +  "overload the " + op + "= operator.\n";

  return str;
}

double AtomicParameter::getValue(int ival) const
{
  const std::string method = "double AtomicParameter::getValue(int) const";
  if(ival >= (int)values.size() && values.size() != 1)
  {
     antiochError(method,"Trying to get out-of-range value.");
     return 0.;
  }

  return (ival >= (int)values.size())?values[0]:values[ival];
  
}

void AtomicParameter::setValue(double val,int ival)
{
  if((values.size() == 0 && ival != 0) || ival > (int)values.size())
  {
     antiochError("void AtomicParameter::setValue(double,int)",
                   "Trying to assign out-of-range value.");
  }
  if(ival == (int)values.size())
        {values.push_back(val);}
  else
        {values[ival] = val;}
}

void AtomicParameter::addValues(const std::vector<double> &vals)
{
  for(unsigned int i = 0; i < vals.size(); i++)
  {
    addValue(vals[i]);
  }
}

void AtomicParameter::duplicateVal(int ival, int ntimes)
{
  const std::string method("void AtomicParameter::duplicateVal(int, int)");
  if(ival >= (int)values.size())
  {
     antiochError(method,"Trying to access out-of-range value.");
  }
  double val = values[ival];
  values.insert(values.begin() + ival,ntimes,val);
}


int AtomicParameter::testArithmOperation(const std::string &method, const AtomicParameter &rhs, 
                                          const std::string &op,     const std::string &sym) const
{
  if(values.empty())return 0;
//dealing with errors: 
//  either one is constant: sizeVal() == 1
//  or two have same number of values
  int min = ((int)values.size() < rhs.sizeVal())?(int)values.size():rhs.sizeVal();
  int max = ((int)values.size() > rhs.sizeVal())?(int)values.size():rhs.sizeVal();
  if( (min != 1) && (min != max))antiochError(method,computationError(op,sym));

  return max;
}

const bool AtomicParameter::operator==(const AtomicParameter &rhs) const
{
  bool isEqual = ((name == rhs.getName()) && ((int)values.size() == rhs.sizeVal()));
  if(isEqual)
  {
    for(unsigned int i = 0; i < values.size(); i++)
    {
      isEqual = ((isEqual) && (rhs.getValue(i)));
    }
  }
  return isEqual;
}

const bool AtomicParameter::operator!=(const AtomicParameter &rhs) const
{
  return (!(*this == rhs));
}

AtomicParameter &AtomicParameter::operator=(const AtomicParameter &rhs)
{
  if(this == &rhs)return *this;
  if(name.empty())name = rhs.getName();
  values = rhs.getValues();
  return *this;
}

AtomicParameter & AtomicParameter::operator+=(const AtomicParameter &rhs)
{//no name considerations
  const std::string method("AtomicParameter & AtomicParameter::operator+=(const AtomicParameter &rhs)");

  int max = testArithmOperation(method,rhs,"add","+");

  int i,z(0.);
  int *k(&i);
  if(sizeVal() != max)values.resize(max,values[0]);
  if(rhs.sizeVal() != max)k = &z;
  for(i = 0; i < max; i++)
  {
     values[i] += rhs.getValue(*k);
  }
  return *this;
}

AtomicParameter & AtomicParameter::operator*=(const AtomicParameter &rhs)
{
  const std::string method("AtomicParameter & AtomicParameter::operator*=(const AtomicParameter &rhs)");
  int max = testArithmOperation(method,rhs,"multiply","*");
  int i,z(0);
  int *k(&i);
  if(max != (int)values.size())values.resize(max,values[0]);
  if(max > rhs.sizeVal())k = &z;
  for(i = 0; i < max; i++)
  {
    values[i] *= rhs.getValue(*k);
  }
  return *this;
}

AtomicParameter & AtomicParameter::operator/=(const AtomicParameter &rhs)
{
  const std::string method("AtomicParameter & AtomicParameter::operator/=(const AtomicParameter &)");
  int max = testArithmOperation(method,rhs,"divid","/");
  int i,z(0);
  int *k(&i);
  if(max != (int)values.size())values.resize(max,values[0]);
  if(max > rhs.sizeVal())k = &z;
  for(i = 0; i < max; i++)
  {
    values[i] /= rhs.getValue(*k);
  }
  return *this;
}

AtomicParameter AtomicParameter::operator+(const AtomicParameter &rhs) const
{
  return (AtomicParameter(*this) += rhs);
}

AtomicParameter AtomicParameter::operator*(const AtomicParameter &rhs) const
{
  return (AtomicParameter(*this) *= rhs);
}

AtomicParameter & AtomicParameter::operator*=(double r)
{
  for(int i = 0; i < sizeVal(); i++)
  {
     values[i] *= r;
  }

  return *this;
}

AtomicParameter AtomicParameter::operator*(double r) const 
{
  return (AtomicParameter(*this) *= r);
}

AtomicParameter & AtomicParameter::operator/=(double r)
{
  *this *= 1./r;
  return *this;
}

AtomicParameter AtomicParameter::operator/(double r) const
{
  return (AtomicParameter(*this) /= r);
}

AtomicParameter & AtomicParameter::operator+=(double r)
{
  for(int i = 0; i < sizeVal(); i++)
  {
     values[i] += r;
  }

  return *this;
}

AtomicParameter AtomicParameter::operator+(double r) const 
{
  return (AtomicParameter(*this) += r);
}

AtomicParameter & AtomicParameter::operator-=(double r)
{
  r *= -1.;
  *this += r;
  return *this;
}

AtomicParameter AtomicParameter::operator-(double r) const 
{
  return (AtomicParameter(*this) -= r);
}

AtomicParameter  AtomicParameter::operator/ (const AtomicParameter &rhs) const
{
  return (AtomicParameter(*this) /= rhs);
}

AtomicParameter & AtomicParameter::operator-=(const AtomicParameter &rhs)
{
  *this += -rhs;
  return *this;
}

AtomicParameter AtomicParameter::operator-(const AtomicParameter &rhs) const
{
  return (AtomicParameter(*this) -= rhs);
}

AtomicParameter AtomicParameter::operator-() const
{
  AtomicParameter out(name);
  for(int i = 0; i < sizeVal(); i++)
  {
     out.addValue(-values[i]);
  }
  return out;
}

}
