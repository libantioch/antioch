//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_EUCKEN_THERMAL_CONDUCTIVITY_UTILS_DECL_H
#define ANTIOCH_EUCKEN_THERMAL_CONDUCTIVITY_UTILS_DECL_H

// Antioch
#include "antioch/physics_metaprogramming_decl.h"

namespace Antioch
{

  // tag
  struct eucken_thermal_conductivity_tag{};

  // getting the tag
  template <typename ThermoEvaluator, typename CoeffType>
  struct physical_tag<EuckenThermalConductivity<ThermoEvaluator, CoeffType> >;

  // not a physical set
   template<typename Model, typename StateType, typename VectorStateType>
   void  physical_set_operator_thermal_conductivity(const Model & set, unsigned int s, const StateType & mu, const StateType & /*dss*/, const StateType & T, const StateType & /*rho*/, 
                                                    StateType & k, eucken_thermal_conductivity_tag);
}

#endif
