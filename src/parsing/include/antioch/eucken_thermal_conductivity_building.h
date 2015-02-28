//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014 Paul T. Bauman, Benjamin S. Kirk, Sylvain Plessis,
//                    Roy H. Stonger
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
#include "antioch/physical_set.h"

// C++
#include <iostream>
#include <vector>


namespace Antioch
{
  //forward declaration
  template <typename ThermoEval, typename NumericType>
  class TransportMixture;

  // backward compatibility
  template<class ThermoEvaluator, class Mixture>
  void build_eucken_thermal_conductivity( PhysicalSet<EuckenThermalConductivity<ThermoEvaluator> , Mixture>& k){}

  template<class ThermoEucken,class ThermoTran, class NumericType>
  void build_eucken_thermal_conductivity( PhysicalSet<EuckenThermalConductivity<ThermoEucken> , TransportMixture<ThermoTran,NumericType> >& k);

// ----------------------------------------- //

  template<class ThermoEucken,class ThermoTran, class NumericType>
  void build_eucken_thermal_conductivity( PhysicalSet<EuckenThermalConductivity<ThermoEucken>, TransportMixture<ThermoTran,NumericType> >& k)
  {
       k.add_model(k.mixture().thermo().micro_thermo());
  }

}

#endif
