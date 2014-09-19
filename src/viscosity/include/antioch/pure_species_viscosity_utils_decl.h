//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_PURE_SPECIES_VISCOSITY_UTILS_DECL_H
#define ANTIOCH_PURE_SPECIES_VISCOSITY_UTILS_DECL_H

#include "antioch/physics_metaprogramming_decl.h"

namespace Antioch
{
   // tag
   struct pure_species_viscosity_tag{};

   // getting tag
   template <typename CoeffType, typename Interpolator>
   struct physical_tag<PureSpeciesViscosity<CoeffType, Interpolator> >;

   // physical set boolean
   template<typename CoeffType, typename Interpolator>
   struct is_physical_set<PureSpeciesViscosity<CoeffType, Interpolator> >;

   // we can initialize without the user's help,
   template <typename ModelSet>
   void physical_set_initialize(ModelSet & mod, pure_species_viscosity_tag );

   // species
   template<typename Model, typename StateType>
   ANTIOCH_AUTO(StateType) 
        physical_set_first_operator(const Model & set, unsigned int s, const StateType & T, pure_species_viscosity_tag);

   // mixture
   template<typename Model, typename StateType, typename VectorStateType>
   void physical_set_operator_viscosity(const Model & set, const StateType & T, VectorStateType & mu, sutherland_viscosity_tag);
}

#endif
