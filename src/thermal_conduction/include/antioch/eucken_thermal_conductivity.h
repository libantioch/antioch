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

#include "antioch/antioch_asserts.h"
#include "antioch/species_conductivity_base.h"

namespace Antioch
{
  template<class ThermoEvaluator>
  class EuckenThermalConductivity : public SpeciesConductivityBase<EuckenThermalConductivity<ThermoEvaluator> >
  {
  public:

    EuckenThermalConductivity( const ThermoEvaluator& t)
      : SpeciesConductivityBase<EuckenThermalConductivity<ThermoEvaluator> >(),
      _thermo(t)
    {}

    template <typename CoeffType>
    EuckenThermalConductivity( const ThermoEvaluator& t, const std::vector<CoeffType>& /*coeffs*/)
      : SpeciesConductivityBase<EuckenThermalConductivity<ThermoEvaluator> >(),
      _thermo(t)
    {}

    virtual ~EuckenThermalConductivity() {}

    template <typename StateType>
    ANTIOCH_AUTO(StateType)
    trans( const unsigned int s, const StateType& mu ) const
    ANTIOCH_AUTOFUNC(StateType,
		     typename Antioch::raw_value_type<StateType>::type(2.5) *
		     mu * this->_thermo.cv_trans(s))

    template <typename StateType>
    ANTIOCH_AUTO(StateType)
    rot( const unsigned int s, const StateType& mu ) const
    ANTIOCH_AUTOFUNC(StateType, mu* this->_thermo.cv_rot(s))

    template <typename StateType>
    ANTIOCH_AUTO(StateType)
    vib( const unsigned int s, const StateType& mu, const StateType& Tv ) const
    /* Note: Cander, "High-temperature effects in hypersonic flight",
             Encyclopedia of Aerospace Engineering, 2010 suggests that
             there should be a factor of 1.2 here. */
    ANTIOCH_AUTOFUNC(StateType, mu* this->_thermo.cv_vib(s, Tv))

    template <typename StateType>
    ANTIOCH_AUTO(StateType)
    elec( const unsigned int s, const StateType& mu, const StateType& Te ) const
    ANTIOCH_AUTOFUNC(StateType, mu* this->_thermo.cv_el(s, Te))

    typedef ThermoEvaluator micro_thermo_type;

    //! Friend the base class so we can make the implementation protected
    friend class SpeciesConductivityBase<EuckenThermalConductivity<ThermoEvaluator> >;

  protected:

    const ThermoEvaluator& _thermo;

     // Suppose thermal equilibrium Tv = Te = T:w
    template <typename StateType>
    ANTIOCH_AUTO(StateType)
      op_no_diff_impl( const unsigned int s, const StateType& mu, const StateType & T ) const
    ANTIOCH_AUTOFUNC(StateType, trans(s,mu) + rot(s,mu) + vib(s,mu,T) + elec(s,mu,T))


    template <typename StateType>
    const ANTIOCH_AUTO(StateType)
      op_with_diff_impl(unsigned int s, const StateType& mu_s, const StateType & T, const StateType & rho_s, const StateType & Dss) const;

  };

  template<class ThermoEvaluator>
  template <typename StateType>
  const ANTIOCH_AUTO(StateType)
    EuckenThermalConductivity<ThermoEvaluator>::op_with_diff_impl(unsigned int s, const StateType& mu_s, const StateType & /*T*/, const StateType & /*rho_s*/, const StateType & /*Dss*/) const
  {
    antioch_error();

    /*The following is dummy*/
    return trans(s,mu_s);
  }

} // end namespace Antioch


#endif // ANTIOCH_EUCKEN_THERMAL_CONDUCTIVITY_H
