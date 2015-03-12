//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-
//--------------------------------------------------------------------------

#ifndef ANTIOCH_PHYSICS_UTILS_H
#define ANTIOCH_PHYSICS_UTILS_H

// Antioch
#include "antioch/physical_set.h"
#include "antioch/transport_mixture.h"

//C++

namespace Antioch
{

  /*! a free fonction to extrapolate to high temperatures
      
     Only the Stockmayer potential imposes a maximal temperature,
     which means only the kinetics theory on transport is
     limited by it.
  */

  // We are explicit on the mixture
  template <typename Physics, typename ThermoEvaluator, typename CoeffType, typename StateType>
  void extrapolate_to_high_temperatures(PhysicalSet<Physics,TransportMixture<ThermoEvaluator,CoeffType> > & physics, const StateType & T_max);


// ---------------------------


  template <typename Physics, typename ThermoEvaluator, typename CoeffType, typename StateType>
  inline
  void extrapolate_to_high_temperatures(PhysicalSet<Physics,TransportMixture<ThermoEvaluator,CoeffType> > & physics, const StateType & T_max)
  {
     extrapolate_T(physics.set(),T_max,typename physical_tag<Physics>::temperature_limitation_type() );
  }

}//end namespace Antioch


#endif //ANTIOCH_PHYSICS_UTILS_H
