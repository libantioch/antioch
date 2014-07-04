//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014 Paul T. Bauman, Sylvain Plessis, Roy H. Stonger
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
#ifndef ANTIOCH_IN_SI_H
#define ANTIOCH_IN_SI_H

//Antioch
#include "antioch/antioch_asserts.h"

//C++
#include <iostream>

namespace Antioch{

/*! \class InSI
 * \brief Seven integers to characterize the power vector
 *
 * A unit is caracterize by the power associated to
 * every base unit, this class stores those power.
 */
class InSI{
    public:

/*! \brief Building constructor, fully descriptive, with zeros as default values to
 * simplify coder's interface*/
      InSI(int i0=0,int i1=0, int i2=0, int i3=0, int i4=0, int i5=0, int i6=0, int i7=0):
        m(i0),kg(i1),s(i2),A(i3),K(i4),mol(i5),cd(i6),rad(i7){}

/*! \brief << operator, to format the power vector*/
      friend std::ostream &operator<< (std::ostream &out, const InSI & rhs)
        {
          out << "["
              << rhs.get_m()   << "(m),"
              << rhs.get_kg()  << "(kg),"
              << rhs.get_s()   << "(s),"
              << rhs.get_A()   << "(A),"
              << rhs.get_K()   << "(K),"
              << rhs.get_mol() << "(mol),"
              << rhs.get_cd()  << "(cd),"
              << rhs.get_rad() << "(rad)]";
          return out;
        }

/*! \brief Bool equalize operator, true if all the powers are equal*/
      bool operator== (const InSI & rhs) const;
/*! \brief Not bool const operator==(const InSI&)*/
      bool operator!= (const InSI & rhs) const;
/*! \brief Assignement operator, equalize all the powers*/
      InSI &     operator=  (const InSI & rhs);
/*! \brief Adding operator, add all the powers*/
      InSI &     operator+= (const InSI & rhs);
/*! \brief Substracting operator, substract all the powers*/
      InSI &     operator-= (const InSI & rhs);
/*! \brief Multiplying operator, multiply all the powers*/
      InSI &     operator*= (int rhs);
/*! \brief Dividing operator.
 *
 * Dividing a power needs first to check that the
 * division is possible, i.e. if the power is not
 * null, it must be a multiple of the divider. If not
 * it sends back an error and stops. This division
 * exist for root purposes.
 */
      InSI &     operator/= (int rhs);
/*! \brief Adding operator, add all the powers*/
      InSI       operator+  (const InSI & rhs) const;
/*! \brief Substracting operator, substract all the powers*/
      InSI       operator-  (const InSI & rhs) const;
/*! \brief Multiplying operator, multiply all the powers*/
      InSI       operator*  (int rhs) const;
/*! \brief Dividing operator, see InSi & operator/=(int) for details*/
      InSI       operator/  (int rhs) const;

/*! \brief Set all the powers to zeros*/
      void clear();

/*! \brief meter power getter*/
      int get_m()   const {return m;}
/*! \brief kilogramme power getter*/
      int get_kg()  const {return kg;}
/*! \brief second power getter*/
      int get_s()   const {return s;}
/*! \brief ampere power getter*/
      int get_A()   const {return A;}
/*! \brief kelvin power getter*/
      int get_K()   const {return K;}
/*! \brief mol power getter*/
      int get_mol() const {return mol;}
/*! \brief candela power getter*/
      int get_cd()  const {return cd;}
/*! \brief radian power getter*/
      int get_rad() const {return rad;}

/*!\brief Check if empty (all values to zero)*/
      bool empty() const;

    private:
      int m,kg,s,A,K,mol,cd,rad;
};

inline
bool InSI::operator== (const InSI & rhs)const
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

inline
bool InSI::operator!= (const InSI & rhs)const
{
  return (!(*this == rhs));
}

inline
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

inline
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

inline
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

inline
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

inline
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

inline
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

inline
InSI & InSI::operator/= (int rhs)
{
  if(m  %rhs != 0)antioch_unit_error("Cannot have non integer power (m).");
  if(kg %rhs != 0)antioch_unit_error("Cannot have non integer power (kg).");
  if(s  %rhs != 0)antioch_unit_error("Cannot have non integer power (s).");
  if(A  %rhs != 0)antioch_unit_error("Cannot have non integer power (A).");
  if(K  %rhs != 0)antioch_unit_error("Cannot have non integer power (K).");
  if(mol%rhs != 0)antioch_unit_error("Cannot have non integer power (mol).");
  if(cd %rhs != 0)antioch_unit_error("Cannot have non integer power (cd).");
  if(rad%rhs != 0)antioch_unit_error("Cannot have non integer power (rad).");
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

inline
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

inline
InSI InSI::operator/ (int rhs)const
{
  return (InSI(*this) /= rhs);
}

inline
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

inline
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
} //Antioch namespace

#endif
