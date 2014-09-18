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
   struct physical_tag_type<ConstantLewisDiffusivity<CoeffType> >;

   template<typename Model, typename StateType, typename ThermoEvaluator>
   void physical_set_operator_diffusion_comes_last(const Model & set, const StateType & cp, const StateType & k, StateType & /*ds*/, constant_lewis_diffusivity_tag);
}

#endif
