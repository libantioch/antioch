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

#ifndef ANTIOCH_SUTHERLAND_VISCOSITY_UTILS_H
#define ANTIOCH_SUTHERLAND_VISCOSITY_UTILS_H

#include "antioch/antioch_asserts.h"
#include "antioch/sutherland_viscosity_utils_decl.h"

namespace Antioch
{
// Sutherland
   // getting tag
   template <typename CoeffType>
   struct physical_tag<SutherlandViscosity<CoeffType> >:
        public physical_tag_base<SutherlandViscosity<CoeffType> >
   {
      typedef sutherland_viscosity_tag type;
        // for operators, diffusion is special, see comment below
     typedef sutherland_viscosity_tag viscosity_type;
   };

   // physical set boolean
   template<typename CoeffType>
   struct is_physical_set<SutherlandViscosity<CoeffType> >
   {
      static const bool value = true;
   };

   template<typename Model, typename StateType>
   void physical_set_operator_viscosity(const Model & set, unsigned int s, const StateType & T, StateType & mu, sutherland_viscosity_tag)
   {
        mu = (*set[s])(T);
   }

   template<typename Model, typename StateType, typename VectorStateType>
   void physical_set_operator_viscosity(const Model & set, const StateType & T, VectorStateType & mu, sutherland_viscosity_tag)
   {
      antioch_assert_equal_to(mu.size(), set.size());

      for(unsigned int s = 0; s < mu.size(); s++)
      {
          mu[s] = (*set[s])(T);
      }
   }
}


#endif
