//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_BLOTTNER_VISCOSITY_UTILS_DECL_H
#define ANTIOCH_BLOTTNER_VISCOSITY_UTILS_DECL_H

#include "antioch/physics_metaprogramming_decl.h"

namespace Antioch
{
// metaprogramming for physics
   // tag
   struct blottner_viscosity_tag{};

   // getting tag
   template <typename CoeffType>
   struct physical_tag_type<BlottnerViscosity<CoeffType> >;

   // physical set boolean
   template<typename CoeffType>
   struct is_physical_set<BlottnerViscosity<CoeffType> >;

   template<typename Model, typename StateType, typename VectorStateType>
   ANTIOCH_AUTO(StateType) 
        physical_set_operator_viscosity(const Model & set, unsigned int s, const StateType & T, blottner_viscosity_tag);
}

#endif
