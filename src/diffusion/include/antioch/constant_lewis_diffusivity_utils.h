//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014 Paul T. Bauman, Benjamin S. Kirk, Sylvain Plessis,
//                    Roy H. Stonger
// Copyright (C) 2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
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
