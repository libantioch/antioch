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

#ifndef ANTIOCH_EUCKEN_THERMAL_CONDUCTIVITY_UTILS_DECL_H
#define ANTIOCH_EUCKEN_THERMAL_CONDUCTIVITY_UTILS_DECL_H

// Antioch
#include "antioch/physics_metaprogramming_decl.h"

namespace Antioch
{

  // tag
  struct eucken_thermal_conductivity_tag{};

  // getting the tag
  template <typename ThermoEvaluator>
  struct physical_tag<EuckenThermalConductivity<ThermoEvaluator> >;

  // we can initialize on our own
   template <typename Model>
   void physical_set_initialize(Model & mod, eucken_thermal_conductivity_tag);

   //species
  template<typename Model, typename StateType>
  void physical_set_operator_thermal_conductivity(const Model & set, unsigned int s, const StateType & mu, const StateType & /*dss*/, const StateType & T, const StateType & /*rho*/, 
                                                   StateType &k, eucken_thermal_conductivity_tag);

  // mixture
  template<typename Model, typename StateType, typename VectorStateType>
  void  physical_set_operator_thermal_conductivity(const Model & set, const VectorStateType & mu, const VectorStateType & /*dss*/, const StateType & T, const StateType & /*rho*/, 
                                                   VectorStateType & k, eucken_thermal_conductivity_tag);
}

#endif
