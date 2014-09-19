//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_CONSTANT_LEWIS_DIFFUSIVITY_UTILS_H
#define ANTIOCH_CONSTANT_LEWIS_DIFFUSIVITY_UTILS_H

#include "antioch/antioch_asserts.h"
#include "antioch/constant_lewis_diffusivity_utils_decl.h"

namespace Antioch
{

   template <typename CoeffType>
   struct physical_tag<ConstantLewisDiffusivity<CoeffType> >:
        public physical_tag_base<ConstantLewisDiffusivity<CoeffType> >
   {
      typedef constant_lewis_diffusivity_tag type;
        // for operators
     typedef constant_lewis_diffusivity_tag diffusion_mixture_type;
   };


   template<typename Model, typename StateType>
   void physical_set_operator_diffusion(const Model & set, const StateType & rho, const StateType & cp, const StateType & k, StateType & ds, constant_lewis_diffusivity_tag)
   {
       antioch_assert(set);
       ds = (*set)(rho,cp,k);
   }

   template<typename Model, typename StateType, typename VectorStateType>
   void physical_set_operator_diffusion(const Model & set, const StateType & rho, const VectorStateType & cp, const VectorStateType & k, VectorStateType & ds, constant_lewis_diffusivity_tag)
   {
       antioch_assert(set);
       antioch_assert_equal_to(k.size(),ds.size());
       antioch_assert(!k.empty());

       for(unsigned int s = 0; s < ds.size(); ++s)
       {
         ds[s] = (*set)(rho,cp[s],k[s]);
       }
   }
}

#endif
