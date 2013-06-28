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
#include "antioch/InSI.hpp"
#include "antioch/Error.hpp"

namespace Antioch{
bool const InSI::operator== (const InSI & rhs)const
{
  return (
           m   == rhs.get_m()   &&
           kg  == rhs.get_kg()  &&
           s   == rhs.get_s()   &&
           A   == rhs.get_A()   &&
           K   == rhs.get_K()   &&
           mol == rhs.get_mol() &&
           cd  == rhs.get_cd()  &&
           rad == rhs.get_rad()
         );
}

bool const InSI::operator!= (const InSI & rhs)const
{
  return (!(*this == rhs));
}

InSI & InSI::operator= (const InSI & rhs)
{
  if(this == &rhs){return *this;}
  m   = rhs.get_m();
  kg  = rhs.get_kg();
  s   = rhs.get_s();
  A   = rhs.get_A();
  K   = rhs.get_K();
  mol = rhs.get_mol();
  cd  = rhs.get_cd();
  rad = rhs.get_rad();
  return *this;
}

InSI InSI::operator+ (const InSI & rhs)const
{
  return InSI(
              m   + rhs.get_m()  ,
              kg  + rhs.get_kg() ,
              s   + rhs.get_s()  ,
              A   + rhs.get_A()  ,
              K   + rhs.get_K()  ,
              mol + rhs.get_mol(),
              cd  + rhs.get_cd() ,
              rad + rhs.get_rad()
             );
}

InSI InSI::operator- (const InSI & rhs)const
{
   return InSI(
                m   - rhs.get_m()  ,
                kg  - rhs.get_kg() ,
                s   - rhs.get_s()  ,
                A   - rhs.get_A()  ,
                K   - rhs.get_K()  ,
                mol - rhs.get_mol(),
                cd  - rhs.get_cd() ,
                rad - rhs.get_rad()
               );
}

InSI & InSI::operator+= (const InSI & rhs)
{
   m   += rhs.get_m();
   kg  += rhs.get_kg();
   s   += rhs.get_s();
   A   += rhs.get_A();
   K   += rhs.get_K();
   mol += rhs.get_mol();
   cd  += rhs.get_cd();
   rad += rhs.get_rad();
   return *this;
}

InSI & InSI::operator-= (const InSI & rhs)
{
  m   -= rhs.get_m();
  kg  -= rhs.get_kg();
  s   -= rhs.get_s();
  A   -= rhs.get_A();
  K   -= rhs.get_K();
  mol -= rhs.get_mol();
  cd  -= rhs.get_cd();
  rad -= rhs.get_rad();
  return *this;
}

InSI & InSI::operator*= (int rhs)
{
  m   *= rhs;
  kg  *= rhs;
  s   *= rhs;
  A   *= rhs;
  K   *= rhs;
  mol *= rhs;
  cd  *= rhs;
  rad *= rhs;
  return *this;
}

InSI & InSI::operator/= (int rhs)
{
  if(m  %rhs != 0)antiochError("InSI &InSI::operator/=(int)","Cannot have non integer power (m).");
  if(kg %rhs != 0)antiochError("InSI &InSI::operator/=(int)","Cannot have non integer power (kg).");
  if(s  %rhs != 0)antiochError("InSI &InSI::operator/=(int)","Cannot have non integer power (s).");
  if(A  %rhs != 0)antiochError("InSI &InSI::operator/=(int)","Cannot have non integer power (A).");
  if(K  %rhs != 0)antiochError("InSI &InSI::operator/=(int)","Cannot have non integer power (K).");
  if(mol%rhs != 0)antiochError("InSI &InSI::operator/=(int)","Cannot have non integer power (mol).");
  if(cd %rhs != 0)antiochError("InSI &InSI::operator/=(int)","Cannot have non integer power (cd).");
  if(rad%rhs != 0)antiochError("InSI &InSI::operator/=(int)","Cannot have non integer power (rad).");
  m   /= rhs;
  kg  /= rhs;
  s   /= rhs;
  A   /= rhs;
  K   /= rhs;
  mol /= rhs;
  cd  /= rhs;
  rad /= rhs;

  return *this;
}

InSI InSI::operator* (int rhs)const
{
  return InSI(
               m * rhs,
               kg * rhs,
               s * rhs,
               A * rhs,
               K * rhs,
               mol * rhs,
               cd * rhs,
               rad * rhs
              );
}

InSI InSI::operator/ (int rhs)const
{
  return (InSI(*this) /= rhs);
}



void InSI::clear()
{
  m   = 0;
  kg  = 0;
  s   = 0;
  A   = 0;
  K   = 0;
  mol = 0;
  cd  = 0;
  rad = 0;
}

bool InSI::empty() const
{
  return (m   == 0 &&
          kg  == 0 &&
          s   == 0 &&
          A   == 0 &&
          K   == 0 &&
          mol == 0 &&
          cd  == 0 &&
          rad == 0);
}
}
