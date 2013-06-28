//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
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
