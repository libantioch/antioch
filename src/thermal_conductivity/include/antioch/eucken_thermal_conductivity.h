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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef ANTIOCH_EUCKEN_THERMAL_CONDUCTIVITY_H
#define ANTIOCH_EUCKEN_THERMAL_CONDUCTIVITY_H

namespace Antioch
{
  template<class ThermoEvaluator>
  class EuckenThermalConductivity
  {
  protected:

    const ThermoEvaluator& _thermo;

  public:

    EuckenThermalConductivity( const ThermoEvaluator& t) : _thermo(t) {}

    ~EuckenThermalConductivity() {}

     // As seen coded in previous WilkeEvaluator, frankly, all Ts the same??
    template <typename StateType>
    ANTIOCH_AUTO(StateType)
    operator()( const unsigned int s, const StateType& mu, const StateType & T ) const
    ANTIOCH_AUTOFUNC(StateType,
		     mu * (typename Antioch::raw_value_type<StateType>::type(2.5) * _thermo.cv_trans(s) +
                            _thermo.cv_rot(s) + _thermo.cv_vib(s,T) + _thermo.cv_el(s,T)) )


    template <typename StateType>
    ANTIOCH_AUTO(StateType)
    trans( const unsigned int s, const StateType& mu ) const
    ANTIOCH_AUTOFUNC(StateType,
		     typename Antioch::raw_value_type<StateType>::type(2.5) *
		     mu * _thermo.cv_trans(s))

    template <typename StateType>
    ANTIOCH_AUTO(StateType)
    rot( const unsigned int s, const StateType& mu ) const
    ANTIOCH_AUTOFUNC(StateType, mu*_thermo.cv_rot(s))

    template <typename StateType>
    ANTIOCH_AUTO(StateType)
    vib( const unsigned int s, const StateType& mu, const StateType& Tv ) const
    /* Note: Cander, "High-temperature effects in hypersonic flight",
             Encyclopedia of Aerospace Engineering, 2010 suggests that
             there should be a factor of 1.2 here. */
    ANTIOCH_AUTOFUNC(StateType, mu*_thermo.cv_vib(s, Tv))

    template <typename StateType>
    ANTIOCH_AUTO(StateType)
    elec( const unsigned int s, const StateType& mu, const StateType& Te ) const
    ANTIOCH_AUTOFUNC(StateType, mu*_thermo.cv_el(s, Te))

  };
} // end namespace Antioch

#include "antioch/eucken_thermal_conductivity_utils_decl.h"

#endif // ANTIOCH_EUCKEN_THERMAL_CONDUCTIVITY_H
