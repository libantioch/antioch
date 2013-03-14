//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// Antioch - A Gas Dynamics Thermochemistry Library
//
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
  public:

    EuckenThermalConductivity( const ThermoEvaluator& thermo );
    ~EuckenThermalConductivity();

    template <typename StateType>
    StateType k_trans( const StateType mu ) const;

    template <typename StateType>
    StateType k_rot( const StateType mu ) const;

    template <typename StateType>
    StateType k_vib( const StateType mu ) const;

    template <typename StateType>
    StateType k_elec( const StateType mu ) const;

  protected:

    const ThermoEvaluator& _thermo;

  };

  EuckenThermalConductivity::EuckenThermalConductivity( const ThermoEvaluator& thermo )
    : _thermo(thermo)
  {
    return;
  }

  EuckenThermalConductivity::~EuckenThermalConductivity()
  {
    return;
  }

  template <typename StateType>
  StateType EuckenThermalConductivity::k_trans( const unsigned int s, const StateType mu ) const
  {
    return 2.5*mu*_thermo.cv_trans(s);
  }

  template <typename StateType>
  StateType EuckenThermalConductivity::k_rot( const unsigned int s, const StateType mu ) const
  {
    return mu*_thermo.cv_rot(s);
  }

  template <typename StateType>
  StateType EuckenThermalConductivity::k_vib( const unsigned int s, const StateType mu ) const
  {
    /* Note: Cander, "High-temperature effects in hypersonic flight",
             Encyclopedia of Aerospace Engineering, 2010 suggests that
             there should be a factor of 1.2 here. */
    return mu*_thermo.cv_vib(s);
  }

  template <typename StateType>
  StateType EuckenThermalConductivity::k_elec( const unsigned int s, const StateType mu ) const
  {
    return mu*_thermo.cv_elec(s);
  }

} // end namespace Antioch

#endif // ANTIOCH_EUCKEN_THERMAL_CONDUCTIVITY_H
