//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
#include "antioch/Converter.hpp"

namespace Antioch{
Converter & Converter::operator=  (const Converter &rhs)
{
  if(this == &rhs){return *this;}
  a = rhs.geta();
  b = rhs.getb();
  return *this;
}

Converter & Converter::operator*= (double coef)
{
  a *= coef;
  return *this;
}

Converter & Converter::operator/= (double coef)
{
  a /= coef;
  return *this;
}

Converter & Converter::operator*= (const Converter &rhs)
{
  a *= rhs.geta();
  b  = (b+rhs.getb())/rhs.geta();
  return *this;
}

Converter & Converter::operator/= (const Converter &rhs)
{
  a /= rhs.geta();
  b  = (b-rhs.getb())/rhs.geta();
  return *this;
}

Converter & Converter::operator+= (const Converter &rhs)
{
  a *= rhs.geta();
  b += rhs.getb();
  return *this;
}

Converter Converter::operator*  (double coef) const
{
  return Converter(a*coef,b);
}

Converter Converter::operator/  (double coef) const
{
  return Converter(a/coef,b);
}

Converter Converter::operator/  (const Converter &rhs) const
{
  return Converter(a/rhs.geta(),(b-rhs.getb())/rhs.geta());
}

Converter Converter::operator*  (const Converter &rhs) const
{
  return Converter(a * rhs.geta() , (b + rhs.getb())/rhs.geta());
}
}
