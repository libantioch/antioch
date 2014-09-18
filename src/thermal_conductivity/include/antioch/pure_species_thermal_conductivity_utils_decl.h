//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_PURE_SPECIES_THERMAL_CONDUCTIVITY_UTILS_DECL_H
#define ANTIOCH_PURE_SPECIES_THERMAL_CONDUCTIVITY_UTILS_DECL_H

// Antioch
#include "antioch/physics_metaprogramming_decl.h"

namespace Antioch
{
   // tag
   struct pure_species_thermal_conductivity_tag{};

   // getting tag
   template <typename CoeffType>
   struct physical_tag_type<PureSpeciesThermalConductivity<CoeffType> >;

   // physical set boolean
   template<typename ThermoEvaluator, typename CoeffType>
   struct is_physical_set<PureSpeciesThermalConductivity<ThermoEvaluator,CoeffType> >;

  template <typename ThermoEvaluator, typename CoeffType>
  struct Initializer<PureSpeciesThermalConductivity<ThermoEvaluator,CoeffType> >;

   // operator
   template<typename Model, typename StateType, typename VectorStateType>
   void physical_set_operator_thermal_conductivity(const Model & set, unsigned int s, const StateType & mu, const StateType & dss, const StateType & T, const StateType & rho, 
                                    StateType & k, pure_species_thermal_conductivity_tag);
   //initialize me
  template <typename ThermoEvaluator, typename CoeffType>
  struct Initializer<PureSpeciesThermalConductivity<ThermoEvaluator,CoeffType> >;

}

#endif