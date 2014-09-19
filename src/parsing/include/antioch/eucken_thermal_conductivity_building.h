//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
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
  template<class ThermoEvaluator, class NumericType>
  void build_eucken_thermal_conductivity( PhysicalSet<EuckenThermalConductivity<ThermoEvaluator> , TransportMixture<ThermoEvaluator,NumericType> >& k);

// ----------------------------------------- //

  template<class ThermoEvaluator, class NumericType>
  void build_eucken_thermal_conductivity( PhysicalSet<EuckenThermalConductivity<ThermoEvaluator>, TransportMixture<ThermoEvaluator,NumericType> >& k)
  {
       k.add_model(k.mixture().thermo());
  }

}

#endif
