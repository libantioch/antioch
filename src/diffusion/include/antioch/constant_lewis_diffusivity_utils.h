//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_CONSTANT_LEWIS_DIFFUSIVITY_UTILS_H
#define ANTIOCH_CONSTANT_LEWIS_DIFFUSIVITY_UTILS_H

#include "antioch/constant_lewis_diffusivity_utils_decl.h"

namespace Antioch
{

   template <typename CoeffType>
   struct physical_tag_type<ConstantLewisDiffusivity<CoeffType> >
   {
      typedef constant_lewis_diffusivity_tag type;
   };


   template<typename Model, typename StateType, typename ThermoEvaluator>
   void physical_set_operator_diffusion_comes_last(const Model & set, const StateType & cp, const StateType & k, StateType & /*ds*/, constant_lewis_diffusivity_tag)
   {
       ds = (*set)(rho,T,k);
   }
}

#endif
