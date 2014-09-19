//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_EUCKEN_THERMAL_CONDUCTIVITY_UTILS_H
#define ANTIOCH_EUCKEN_THERMAL_CONDUCTIVITY_UTILS_H

#include "antioch/antioch_asserts.h"
#include "antioch/eucken_thermal_conductivity_utils_decl.h"
#include "antioch/eucken_thermal_conductivity_building.h"

namespace Antioch
{
/// Eucken

  // getting the tag
  template <typename ThermoEvaluator>
  struct physical_tag<EuckenThermalConductivity<ThermoEvaluator> >
  {
      typedef EuckenThermalConductivity<ThermoEvaluator> Model;
      typedef eucken_thermal_conductivity_tag type;
        // kind of set tag
     typedef typename if_else_type<is_physical_set<Model>::value,
                                     default_physical_set_tag,
                                     default_physical_tag
                                   >::type set_type;
        // some models require specific initialization
        // (typically automatic initialization)
     typedef eucken_thermal_conductivity_tag init_type;
        // some models require specific initialization
        // but not specific deletion
     typedef typename if_else_type<is_physical_set<Model>::value,
                                     default_physical_set_tag,
                                     default_physical_tag
                                   >::type del_type;
        // for operators
     typedef default_physical_tag            viscosity_type;
     typedef default_physical_tag            diffusion_species_type;
     typedef default_physical_tag            diffusion_mixture_type;
     typedef eucken_thermal_conductivity_tag thermal_conductivity_type;
  };


   template <typename Model>
   void physical_set_initialize(Model & mod, eucken_thermal_conductivity_tag)
   {
      mod.set() = NULL;
      build_eucken_thermal_conductivity( mod);
   }


  template<typename Model, typename StateType>
  void physical_set_operator_thermal_conductivity(const Model & set, unsigned int s, const StateType & mu, const StateType & /*dss*/, const StateType & T, const StateType & /*rho*/, 
                                                   StateType & k, eucken_thermal_conductivity_tag)
  {
     antioch_assert(set);
     k = (*set)(s,mu,T);
  }

  template<typename Model, typename StateType, typename VectorStateType>
  void  physical_set_operator_thermal_conductivity(const Model & set, const VectorStateType & mu, const VectorStateType & /*dss*/, const StateType & T, const StateType & /*rho*/, 
                                                   VectorStateType & k, eucken_thermal_conductivity_tag)
  {
     antioch_assert_equal_to(k.size(),set.size());
     for(unsigned int s = 0; s < k.size(); ++s)
     {
        k[s] = (*set[s])(s,mu[s],T);
     }
  }

}

#endif
