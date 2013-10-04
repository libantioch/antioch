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
#ifndef _UNITS_DEFINITIONS_
#define _UNITS_DEFINITIONS_

/*!\file unit_defs.hpp
 * \brief Provided units and prefixes
 *
 * This file contains the known units to be
 * used for unit management. Some physical constant
 * (double form) are also defined. The known units
 * begins with the SI basis:
 * \f[
 *      [\mathrm{meter,kilogram,second,ampere,kelvin,mole,candela,radian}]
 * \f]
 * 
 * The definition of a unit uses the fully-descriptive constructor
 * Units::UnitBase<T>(string sym,string na,double conva,double convb,int mi,int kgi=0, int si=0, int Ai=0, int Ki=0, int moli=0, int cdi=0, int radi=0);
 * Any unit is thus characterized as follow:
 *
 * UnitBase<T>("mmHg","millimeter of mercury",133.322387415,0.,-1,1,-2).
 *
 * The prefixes are then defined, in the form
 *
 * SIPrefixes<T>("mu",1e-6).
 *
 */

//Antioch
#include "antioch/unit_store.h"

//C++

namespace Antioch{

namespace UnitBaseStorage{

extern const UnitBaseConstant::UnitBaseStore<long double> storage_unit;
extern const UnitBaseConstant::SIPrefixeStore<long double> storage_prefixe;

inline
const UnitBaseConstant::UnitBaseStore<long double>  known_units()
{
   return storage_unit;
}

inline
const UnitBaseConstant::SIPrefixeStore<long double> known_prefixes()
{
  return storage_prefixe;
}

/*
template <typename T>
inline
void allKnownUnitBase(std::ostream &out = std::cout);

template <typename Statetype>
inline
void allKnownPrefixes(std::ostream &out = std::cout);

template <typename T>
void check_unit_consistency();

template <T>
void allKnownHomogeneousUnitBase(const std::string &target, std::ostream &out = std::cout);
*/
}//end namespace UnitBaseStorage


}//end namespace Antioch
#endif
