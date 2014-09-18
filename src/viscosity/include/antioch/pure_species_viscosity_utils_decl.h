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
   struct physical_tag_type<PureSpeciesViscosity<CoeffType, Interpolator> >;

   // getting other tag
   template <typename CoeffType, typename Interpolator>
   struct physical_set_tag<PureSpeciesViscosity<CoeffType,Interpolator> >;

   // physical set boolean
   template<typename CoeffType, typename Interpolator>
   struct is_physical_set<PureSpeciesViscosity<CoeffType, Interpolator> >;

   // we can initialize without the user's help,
   // so we tag the physical_set_tag:
   // we need to define 
   // physical_set_initialize
   // physical_set_delete
   // physical_set_add
   // physical_set_reset
   // but not the default ones that are not concerned
   template <typename ModelSet>
   void physical_set_initialize(ModelSet & mod, pure_species_viscosity_tag );

   template <typename ModelSet>
   void physical_set_delete(ModelSet & mod, pure_species_viscosity_tag );

   template <typename Model, typename InitType>
   void physical_set_add(unsigned int s, SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, pure_species_viscosity_tag);

   template <typename Model, typename InitType>
   void physical_set_add(unsigned int s, SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, pure_species_viscosity_tag);

   template<typename Model, typename StateType, typename VectorStateType>
   ANTIOCH_AUTO(StateType) 
        physical_set_first_operator(const Model & set, unsigned int s, const StateType & T, pure_species_viscosity_tag);
}

#endif
