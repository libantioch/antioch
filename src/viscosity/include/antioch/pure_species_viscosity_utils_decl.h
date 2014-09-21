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

#ifndef ANTIOCH_PURE_SPECIES_VISCOSITY_UTILS_DECL_H
#define ANTIOCH_PURE_SPECIES_VISCOSITY_UTILS_DECL_H

#include "antioch/physics_metaprogramming_decl.h"

namespace Antioch
{
   // tag
   struct pure_species_viscosity_tag{};

   // getting tag
   template <typename CoeffType, typename Interpolator>
   struct physical_tag<PureSpeciesViscosity<CoeffType, Interpolator> >;

   // physical set boolean
   template<typename CoeffType, typename Interpolator>
   struct is_physical_set<PureSpeciesViscosity<CoeffType, Interpolator> >;

   // we can initialize without the user's help,
   template <typename ModelSet>
   void physical_set_initialize(ModelSet & mod, pure_species_viscosity_tag );

   // species
   template<typename Model, typename StateType>
   ANTIOCH_AUTO(StateType) 
        physical_set_first_operator(const Model & set, unsigned int s, const StateType & T, pure_species_viscosity_tag);

   // mixture
   template<typename Model, typename StateType, typename VectorStateType>
   void physical_set_operator_viscosity(const Model & set, const StateType & T, VectorStateType & mu, pure_species_viscosity_tag);
}

#endif
