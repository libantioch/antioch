//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_PURE_SPECIES_THERMAL_CONDUCTIVITY_UTILS_H
#define ANTIOCH_PURE_SPECIES_THERMAL_CONDUCTIVITY_UTILS_H

#include "antioch/pure_species_thermal_conductivity_utils_decl.h"

namespace Antioch
{
   // getting tag
   template <typename CoeffType>
   struct physical_tag_type<PureSpeciesThermalConductivity<CoeffType> >
   {
      typedef pure_species_thermal_conductivity_tag type;
   };

   // physical set boolean
   template<typename ThermoEvaluator, typename CoeffType>
   struct is_physical_set<PureSpeciesThermalConductivity<ThermoEvaluator,CoeffType> >
   {
      static const bool value = true;
   };

   // operator
   template<typename Model, typename StateType, typename VectorStateType>
   void physical_set_operator_thermal_conductivity(const Model & set, unsigned int s, const StateType & mu, const StateType & dss, const StateType & T, const StateType & rho, 
                                    StateType & k, pure_species_thermal_conductivity_tag)
   {
        k = set[s](mu,T,rho,dss);
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

}

#endif
