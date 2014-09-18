//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_EUCKEN_THERMAL_CONDUCTIVITY_UTILS_H
#define ANTIOCH_EUCKEN_THERMAL_CONDUCTIVITY_UTILS_H

#include "antioch/eucken_thermal_conductivity_utils_decl.h"

namespace Antioch
{
/// Eucken

  // getting the tag
  template <typename ThermoEvaluator, typename CoeffType>
  struct physical_tag<EuckenThermalConductivity<ThermoEvaluator, CoeffType> >
  {
      typedef eucken_thermal_conductivity_tag type;
  };

  template<typename Model, typename StateType, typename VectorStateType>
  void  physical_set_operator_thermal_conductivity(const Model & set, unsigned int s, const StateType & mu, const StateType & /*dss*/, const StateType & T, const StateType & /*rho*/, 
                                                   StateType & k, eucken_thermal_conductivity_tag)
  {
     k = (*set)(s,mu,T);
  }

}

#endif
