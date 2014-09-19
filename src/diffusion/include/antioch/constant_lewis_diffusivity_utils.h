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
   struct physical_tag<ConstantLewisDiffusivity<CoeffType> >
   {
      typedef ConstantLewisDiffusivity<CoeffType> Model;
      typedef constant_lewis_diffusivity_tag type;
        // kind of set tag
     typedef typename if_else_type<is_physical_set<Model>::value,
                                     default_physical_set_tag,
                                     default_physical_tag
                                   >::type set_type;
        // some models require specific initialization
        // (typically automatic initialization)
     typedef typename if_else_type<is_physical_set<Model>::value,
                                     default_physical_set_tag,
                                     default_physical_tag
                                   >::type init_type;
        // some models require specific initialization
        // but not specific deletion
     typedef typename if_else_type<is_physical_set<Model>::value,
                                     default_physical_set_tag,
                                     default_physical_tag
                                   >::type del_type;
        // for operators
     typedef default_physical_tag           viscosity_type;
     typedef default_physical_tag           diffusion_species_type;
     typedef constant_lewis_diffusivity_tag diffusion_mixture_type;
     typedef default_physical_tag           thermal_conductivity_type;
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
