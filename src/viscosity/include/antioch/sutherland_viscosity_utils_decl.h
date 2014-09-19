//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_SUTHERLAND_VISCOSITY_UTILS_DECL_H
#define ANTIOCH_SUTHERLAND_VISCOSITY_UTILS_DECL_H

#include "antioch/physics_metaprogramming_decl.h"

namespace Antioch
{
// metaprogramming for physics
   // tag
   struct sutherland_viscosity_tag{};

   // getting tag
   template <typename CoeffType>
   struct physical_tag<SutherlandViscosity<CoeffType> >;

   // physical set boolean
   template<typename CoeffType>
   struct is_physical_set<SutherlandViscosity<CoeffType> >;

   // species
   template<typename Model, typename StateType>
   void physical_set_operator_viscosity(const Model & set, unsigned int s, const StateType & T, StateType & mu, sutherland_viscosity_tag);

   // mixture
   template<typename Model, typename StateType, typename VectorStateType>
   void physical_set_operator_viscosity(const Model & set, const StateType & T, VectorStateType & mu, sutherland_viscosity_tag);
}

#endif
