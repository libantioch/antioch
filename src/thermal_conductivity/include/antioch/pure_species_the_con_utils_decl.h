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

#ifndef ANTIOCH_PURE_SPECIES_THERMAL_CONDUCTIVITY_UTILS_DECL_H
#define ANTIOCH_PURE_SPECIES_THERMAL_CONDUCTIVITY_UTILS_DECL_H

// Antioch
#include "antioch/physics_metaprogramming_decl.h"

namespace Antioch
{
   // tag
   struct pure_species_thermal_conductivity_tag{};

   // getting tag
   template <typename ThermoEvaluator, typename CoeffType>
   struct physical_tag<PureSpeciesThermalConductivity<ThermoEvaluator, CoeffType> >;

   // physical set boolean
   template<typename ThermoEvaluator, typename CoeffType>
   struct is_physical_set<PureSpeciesThermalConductivity<ThermoEvaluator,CoeffType> >;

   //initialize me
   template <typename ThermoEvaluator, typename CoeffType>
   struct Initializer<PureSpeciesThermalConductivity<ThermoEvaluator,CoeffType>, pure_species_thermal_conductivity_tag >;

   // operator
   template<typename Model, typename StateType>
   void physical_set_operator_thermal_conductivity(const Model & set, unsigned int s, const StateType & mu, const StateType & dss, const StateType & T, const StateType & rho, 
                                                   StateType & k, pure_species_thermal_conductivity_tag);

   template<typename Model, typename StateType, typename VectorStateType>
   void physical_set_operator_thermal_conductivity(const Model & set, const VectorStateType & mu, const VectorStateType & dss, const StateType & T, const StateType & rho, 
                                                   VectorStateType & k, pure_species_thermal_conductivity_tag);

   template <typename ModelSet>
   void physical_set_initialize(ModelSet & mod, pure_species_thermal_conductivity_tag);

   // custom add
   template <typename Model, typename InitType>
   void physical_set_add(unsigned int s, typename SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, pure_species_thermal_conductivity_tag);

   // requires to void this
   template <typename Model, typename InitType>
   void physical_set_add(Model & set, const InitType & init, pure_species_thermal_conductivity_tag);

   // custom reset
   template <typename Model, typename InitType>
   void physical_set_reset(unsigned int s, typename SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init,  pure_species_thermal_conductivity_tag);

   // requires to void this
   template <typename Model, typename InitType>
   void physical_set_reset(Model & set, const InitType & init,  pure_species_thermal_conductivity_tag);
}

#endif
