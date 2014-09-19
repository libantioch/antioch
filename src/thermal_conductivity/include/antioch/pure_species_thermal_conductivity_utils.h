//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_PURE_SPECIES_THERMAL_CONDUCTIVITY_UTILS_H
#define ANTIOCH_PURE_SPECIES_THERMAL_CONDUCTIVITY_UTILS_H

#include "antioch/antioch_asserts.h"
#include "antioch/pure_species_thermal_conductivity_utils_decl.h"
#include "antioch/pure_species_thermal_conductivity_building.h"

namespace Antioch
{
   // getting tag
   template <typename CoeffType>
   struct physical_tag<PureSpeciesThermalConductivity<CoeffType> >
   {
      typedef PureSpeciesThermalConductivity<CoeffType> Model;
      typedef pure_species_thermal_conductivity_tag type;
        // kind of set tag
     typedef typename if_else_type<is_physical_set<Model>::value,
                                     default_physical_set_tag,
                                     default_physical_tag
                                   >::type set_type;
        // some models require specific initialization
        // (typically automatic initialization)
     typedef pure_species_thermal_conductivity init_type;
        // some models require specific initialization
        // but not specific deletion
     typedef typename if_else_type<is_physical_set<Model>::value,
                                     default_physical_set_tag,
                                     default_physical_tag
                                   >::type del_type;
        // for operators, diffusion is special, see comment below
     typedef default_physical_tag viscosity_type;
     typedef default_physical_tag diffusion_species_type;
     typedef default_physical_tag diffusion_mixture_type;
     typedef default_physical_tag thermal_conductivity_type;
   };

   // physical set boolean
   template<typename ThermoEvaluator, typename CoeffType>
   struct is_physical_set<PureSpeciesThermalConductivity<ThermoEvaluator,CoeffType> >
   {
      static const bool value = true;
   };

   // operator
   template<typename Model, typename StateType>
   void physical_set_operator_thermal_conductivity(const Model & set, unsigned int s, const StateType & mu, const StateType & dss, const StateType & T, const StateType & rho, 
                                                   StateType & k, pure_species_thermal_conductivity_tag)
   {
        antioch_assert(!set.empty());

        k = (*set[s])(mu,T,rho,dss);
   }

   template<typename Model, typename StateType, typename VectorStateType>
   void physical_set_operator_thermal_conductivity(const Model & set, const VectorStateType & mu, const VectorStateType & dss, const StateType & T, const StateType & rho, 
                                                   VectorStateType & k, pure_species_thermal_conductivity_tag)
   {
      antioch_assert_equal_to(k.size(),set.size());

      for(unsigned int s = 0; s < k.size(); ++s)
      {
          k[s] = (*set[s])(mu[s],T,rho,dss[s]);
      }
   }

   //initialize me
  template <typename ThermoEvaluator, typename CoeffType>
  struct Initializer<PureSpeciesThermalConductivity<ThermoEvaluator,CoeffType> >
  {
     const ThermoEvaluator & t;
     const CoeffType & Z_298K;
     const CoeffType & M;
     const CoeffType & LJ_depth;
  };

   template <typename ModelSet>
   void physical_set_initialize(ModelSet & mod, pure_species_thermal_conductivity_tag)
   {
      mod.set().resize(mod.mixture().n_species(),NULL);
      build_pure_species_thermal_conductivity(mod);
   }

}

#endif
