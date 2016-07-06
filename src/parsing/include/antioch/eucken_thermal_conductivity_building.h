//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014-2016 Paul T. Bauman, Benjamin S. Kirk,
//                         Sylvain Plessis, Roy H. Stonger
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

#ifndef ANTIOCH_EUCKEN_THERMAL_CONDUCTIVITY_BUILDING_H
#define ANTIOCH_EUCKEN_THERMAL_CONDUCTIVITY_BUILDING_H

// Antioch
#include "antioch/eucken_thermal_conductivity.h"
#include "antioch/mixture_conductivity.h"

// C++
#include <iostream>
#include <vector>


namespace Antioch
{
  template<class MicroThermo, class NumericType>
  void build_eucken_thermal_conductivity( MixtureConductivity<EuckenThermalConductivity<MicroThermo>,NumericType>& k,
                                          const MicroThermo& thermo )
  {
    for(unsigned int s = 0; s < k.mixture().n_species(); s++)
       {
         // Eucken doesn't have any coefficients to cache so we provide a dummy vector
         std::vector<NumericType> dummy;
         k.add(s,dummy,thermo);
       }
  }

} // end namespace Antioch

#endif // ANTIOCH_EUCKEN_THERMAL_CONDUCTIVITY_BUILDING_H
