//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_EUCKEN_THERMAL_CONDUCTIVITY_BUILDING_H
#define ANTIOCH_EUCKEN_THERMAL_CONDUCTIVITY_BUILDING_H

// Antioch
#include "antioch/eucken_thermal_conductivity.h"
#include "antioch/tranport_mixture.h"
#include "antioch/physical_set.h"

// C++
#include <iostream>
#include <vector>


namespace Antioch
{
  template<class ThermoEvaluator, class NumericType>
  void build_eucken_thermal_conductivity( PhysicalSet<PureSpeciesThermalConductivity<ThermoEvaluator> , TransportMixture<NumericType> >& k)

// ----------------------------------------- //

  template<calss ThermoEvaluator, class NumericType>
  void build_eucken_thermal_conductivity( PhysicalSet<EuckenThermalConductivity<ThermoEvaluator>, TransportMixture<NumericType> >& k)
  {
       k.set() = new EuckenThermalConductivity(k.mixture().thermo());
  }

}

#endif
