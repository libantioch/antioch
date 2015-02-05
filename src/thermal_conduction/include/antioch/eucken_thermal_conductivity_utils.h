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

#ifndef ANTIOCH_EUCKEN_THERMAL_CONDUCTIVITY_UTILS_H
#define ANTIOCH_EUCKEN_THERMAL_CONDUCTIVITY_UTILS_H

#include "antioch/antioch_asserts.h"
#include "antioch/eucken_thermal_conductivity_utils_decl.h"
#include "antioch/eucken_thermal_conductivity_building.h"

namespace Antioch
{
/// Eucken

  // getting the tag
  template <typename ThermoEvaluator>
  struct physical_tag<EuckenThermalConductivity<ThermoEvaluator> >:
        public physical_tag_base<EuckenThermalConductivity<ThermoEvaluator> >
  {
      typedef eucken_thermal_conductivity_tag type;
        // some models require specific initialization
        // (typically automatic initialization)
     typedef eucken_thermal_conductivity_tag init_type;
        // for operators
     typedef eucken_thermal_conductivity_tag thermal_conductivity_type;
  };


   template <typename Model>
   void physical_set_initialize(Model & mod, eucken_thermal_conductivity_tag)
   {
      mod.set() = NULL;
      build_eucken_thermal_conductivity( mod);
   }


  template<typename Model, typename StateType>
  void physical_set_operator_thermal_conductivity(const Model & set, unsigned int s, const StateType & mu, const StateType & /*dss*/, const StateType & T, const StateType & /*rho*/, 
                                                   StateType & k, eucken_thermal_conductivity_tag)
  {
     antioch_assert(set);
     k = (*set)(s,mu,T);
  }

  template<typename Model, typename StateType, typename VectorStateType>
  void  physical_set_operator_thermal_conductivity(const Model & set, const VectorStateType & mu, const VectorStateType & /*dss*/, const StateType & T, const StateType & /*rho*/, 
                                                   VectorStateType & k, eucken_thermal_conductivity_tag)
  {
     antioch_assert_equal_to(k.size(),set.size());
     for(unsigned int s = 0; s < k.size(); ++s)
     {
        k[s] = (*set[s])(s,mu[s],T);
     }
  }

}

#endif
