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

#ifndef ANTIOCH_PURE_SPECIES_VISCOSITY_UTILS_H
#define ANTIOCH_PURE_SPECIES_VISCOSITY_UTILS_H

#include "antioch/antioch_asserts.h"
#include "antioch/pure_species_viscosity_utils_decl.h"
#include "antioch/pure_species_viscosity_building.h"

namespace Antioch
{
   // getting tag
   template <typename CoeffType, typename Interpolator>
   struct physical_tag<PureSpeciesViscosity<CoeffType, Interpolator> >:
        public physical_tag_base<PureSpeciesViscosity<CoeffType, Interpolator> >
   {
      typedef pure_species_viscosity_tag type;
        // some models require specific initialization
        // (typically automatic initialization)
     typedef pure_species_viscosity_tag init_type;
        // for operators, diffusion is special, see comment below
     typedef pure_species_viscosity_tag viscosity_type;
   };

   // physical set boolean
   template<typename CoeffType, typename Interpolator>
   struct is_physical_set<PureSpeciesViscosity<CoeffType, Interpolator> >
   {
      static const bool value = true;
   };

   // we can initialize without the user's help
   template <typename ModelSet>
   void physical_set_initialize(ModelSet & mod, pure_species_viscosity_tag )
   {
      mod.set().resize(mod.mixture().n_species(),NULL);
      build_pure_species_viscosity(mod);
   }

   template<typename Model, typename StateType>
   void physical_set_operator_viscosity(const Model & set, unsigned int s, const StateType & T, StateType & mu, pure_species_viscosity_tag)
   {
      mu = (*set[s])(T);
   }

   template<typename Model, typename StateType, typename VectorStateType>
   void physical_set_operator_viscosity(const Model & set, const StateType & T, VectorStateType & mu, pure_species_viscosity_tag)
   {
      antioch_assert_equal_to(mu.size(), set.size());

      for(unsigned int s = 0; s < mu.size(); s++)
      {
          mu[s] = (*set[s])(T);
      }
   }

}


#endif
