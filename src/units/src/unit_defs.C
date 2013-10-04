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

//Antioch
#include "antioch/unit_defs.h"

//C++
#include <sstream>
namespace Antioch{

namespace UnitBaseStorage{

const UnitBaseConstant::UnitBaseStore<long double> storage_unit;
const UnitBaseConstant::SIPrefixeStore<long double> storage_prefixe;

/*
inline
const UnitBaseConstant::UnitBaseStore<long double> known_units()
{
   return storage_unit;
}

inline
const UnitBaseConstant::SIPrefixeStore<long double> known_prefixes()
{
  return storage_prefixe;
}

template <typename T>
inline
void allKnownUnitBase(std::ostream &out = std::cout)
{
  out << std::left;
  out << std::setw(15) << "Unit symbol";
  out << std::setw(25) << "Unit name";
  out << std::setw(20) << "Coefficient";   
  out << std::setw(50) << "Projection" 
      << std::endl     << "--" << std::endl;
  for(int it = 0; it < known_units().n_known_units(); it++)
  {
     std::ostringstream symb,name,coef,power;
     symb  << known_units().stored(it).symbol();
     name  << known_units().stored(it).name();
     coef  << known_units().stored(it).converter();
     power << known_units().stored(it).power_array();

     out << std::setw(15) << symb.str();
     out << std::setw(25) << name.str();
     out << std::setw(20) << coef.str();
     out << std::setw(50) << power.str() << std::endl;
  }

  return;
}

template <typename Statetype>
inline
void allKnownPrefixes(std::ostream &out = std::cout)
{
  out << std::left;
  out << std::setw(18) << "Prefixe symbol";
  out << std::setw(18) << "Prefixe name";
  out << std::setw(15) << "Value" 
      << std::endl << "--" << std::endl;
  for(int it = 0; it < known_prefixes().n_known_prefixes(); it++)
  {
     out << std::setw(18) << known_prefixes().stored(it).symbol()
         << std::setw(18) << known_prefixes().stored(it).name()
         << std::setw(15) << known_prefixes().stored(it).value()
         << std::endl;
  }

  return ;
}

template <typename T>
inline
void check_unit_consistency()
{
  for(int i = 0; i < known_units().n_known_units(); i++)
  {
     for(int p = 0; p < known_prefixes().n_known_prefixes(); p++)
     {
        std::string currentU = known_prefixes().stored(p).symbol() + known_prefixes().stored(i).symbol();
        for(int j = 0; j < known_units().n_known_units(); j++)
        {
           for(int q = 0; q < known_prefixes().n_known_prefixes(); q++)
           {
             if(j == i && p == q)continue;
              std::string testU = known_prefixes().stored(q).symbol() + known_units().stored(j).symbol();
              if(currentU == testU)
                antioch_unit_error("The file unit_store.h is not good, check your UnitBaseConstant::UnitBaseStore and UnitBaseConstant::SIPrefixeStore.");
           }
        }
     }  
  }
}

template <T>
inline
void allKnownHomogeneousUnitBase(const std::string &target, std::ostream &out = std::cout)
{
  out << std::left;
  out << "Homogeneous known units to " 
      << target        << std::endl;
  out << std::setw(15) << "Unit symbol";
  out << std::setw(25) << "Unit name";
  out << std::setw(20) << "Coefficient";   
  out << std::setw(50) << "Projection" 
      << std::endl     << "--" << std::endl;
  for(int i = 0; i < ; i++)
  {
     if(!UnitBase::knownUnits[i].isHomogeneous(target))continue;
     std::ostringstream symb,name,coef,power;
     symb  << UnitBase::knownUnits[i].getSymbol();
     name  << UnitBase::knownUnits[i].getName();
     coef  << UnitBase::knownUnits[i].getSIcoef();
     power << UnitBase::knownUnits[i].getPower();

     out << std::setw(15) << symb.str();
     out << std::setw(25) << name.str();
     out << std::setw(20) << coef.str();
     out << std::setw(50) << power.str() << std::endl;
  }

  return out.str();
}
*/
}//end namespace UnitBaseStorage


}//end namespace Antioch
