//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_PURE_SPECIES_THERMAL_CONDUCTIVITY_BUILDING_H
#define ANTIOCH_PURE_SPECIES_THERMAL_CONDUCTIVITY_BUILDING_H

// Antioch
#include "antioch/pure_species_thermal_conductivity.h"
#include "antioch/tranport_mixture.h"
#include "antioch/physical_set.h"

// C++
#include <iostream>
#include <vector>


namespace Antioch
{
  template<class ThermoEvaluator, class NumericType>
  void build_pure_species_thermal_conductivity( PhysicalSet<PureSpeciesThermalConductivity<ThermoEvaluator,NumericType> , TransportMixture<NumericType> >& k)

// ----------------------------------------- //

  template<class NumericType>
  void build_pure_species_thermal_conductivity( PhysicalSet<PureSpeciesThermalConductivity<ThermoEvaluator,NumericType>, TransportMixture<NumericType> >& k)
  {
       for(unsigned int s = 0; s < k.mixture().n_species(); s++)
       {
           Initializer<PureSpeciesThermalConductivity<ThermoEvaluator,NumericType> > init;
           init.t        = k.mixture().thermo();
           init.Z_298K   = k.mixture().transport_species()[s].rotational_relaxation();
           init.M        = k.mixture().transport_species()[s].M();
           init.LJ_depth = k.mixture().transport_species()[s].LJ_depth();
           k.add_model(k.mixture().species_inverse_name_map().at(s),init);
       }
  }

}

#endif
