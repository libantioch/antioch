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
   struct physical_tag<PureSpeciesThermalConductivity<CoeffType> >:
        public physical_tag_base<PureSpeciesThermalConductivity<CoeffType> >
   {
      typedef pure_species_thermal_conductivity_tag type;
        // some models require specific initialization
        // (typically automatic initialization)
     typedef pure_species_thermal_conductivity init_type;
        // for operators, diffusion is special, see comment below
     typedef pure_species_thermal_conductivity thermal_conductivity_type;
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
