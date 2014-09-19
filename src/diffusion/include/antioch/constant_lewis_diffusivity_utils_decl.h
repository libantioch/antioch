//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_CONSTANT_LEWIS_DIFFUSIVITY_UTILS_DECL_H
#define ANTIOCH_CONSTANT_LEWIS_DIFFUSIVITY_UTILS_DECL_H

#include "antioch/physics_metaprogramming_decl.h"

namespace Antioch
{
   // tag
   struct constant_lewis_diffusivity_tag{};

   // getting tag
   template <typename CoeffType>
   struct physical_tag<ConstantLewisDiffusivity<CoeffType> >;

   template<typename Model, typename StateType>
   void physical_set_operator_diffusion(const Model & set, const StateType & rho, const StateType & cp, const StateType & k, StateType & ds, constant_lewis_diffusivity_tag);

   template<typename Model, typename StateType, typename VectorStateType>
   void physical_set_operator_diffusion(const Model & set, const StateType & rho, const VectorStateType & cp, const VectorStateType & k, VectorStateType & ds, constant_lewis_diffusivity_tag);
}

#endif
