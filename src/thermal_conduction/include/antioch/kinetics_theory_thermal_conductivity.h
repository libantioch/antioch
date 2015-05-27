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

#ifndef ANTIOCH_KINETICS_THEORY_THERMAL_CONDUCTIVITY_H
#define ANTIOCH_KINETICS_THEORY_THERMAL_CONDUCTIVITY_H

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/metaprogramming_decl.h"
#include "antioch/rotational_relaxation.h"
#include "antioch/species_conductivity_base.h"

//C++
#include <vector>

namespace Antioch{

  template <typename ThermoEvaluator, typename CoeffType>
  class KineticsTheoryThermalConductivity : public SpeciesConductivityBase<KineticsTheoryThermalConductivity<ThermoEvaluator,CoeffType> >
  {
  public:

        KineticsTheoryThermalConductivity(const ThermoEvaluator & t, CoeffType Z_298K, CoeffType LJ_depth);

        KineticsTheoryThermalConductivity( const ThermoEvaluator& t, const std::vector<CoeffType>& coeffs);

        ~KineticsTheoryThermalConductivity();

        void reset_coeffs(const std::vector<CoeffType> & coeffs);

        template <typename StateType>
        ANTIOCH_AUTO(StateType)
        trans( const unsigned int s, const StateType& mu, const StateType & T, const StateType & rho, const StateType & Dss) const
        ANTIOCH_AUTOFUNC(StateType, mu * this->_thermo.cv_trans(s) * five_over_two * (one - two_over_pi * this->_thermo.cv_rot_over_R(s) / this->_thermo.cv_trans_over_R(s) * this->A(rho * Dss / mu) / this->B(s, T, rho * Dss / mu)))

        template <typename StateType>
        ANTIOCH_AUTO(StateType)
        rot( const unsigned int s, const StateType& mu, const StateType & T, const StateType & rho, const StateType & Dss) const
        ANTIOCH_AUTOFUNC(StateType, rho * Dss * this->_thermo.cv_rot(s) * (one + two_over_pi * this->A(rho * Dss / mu) / this->B(s, T, rho * Dss / mu)))

        template <typename StateType>
        ANTIOCH_AUTO(StateType)
        vib( const unsigned int s, const StateType& mu, const StateType & T, const StateType & rho, const StateType & Dss) const
        ANTIOCH_AUTOFUNC(StateType, rho * Dss * this->_thermo.cv_vib(s,T))

        const ThermoEvaluator& thermo() const
        { return this->_thermo; }

        //! const ref to the rotational relaxation
        const RotationalRelaxation<CoeffType> & rot() const {return _rot;}

    //! Friend the base class so we can make the implementation protected
    friend class SpeciesConductivityBase<KineticsTheoryThermalConductivity<ThermoEvaluator,CoeffType> >;

  protected:

    const ThermoEvaluator& _thermo;

    template <typename StateType>
    const ANTIOCH_AUTO(StateType)
      op_with_diff_impl(unsigned int s, const StateType& mu_s, const StateType & T, const StateType & rho_s, const StateType & Dss) const;

    template <typename StateType>
    ANTIOCH_AUTO(StateType)
      op_no_diff_impl( const unsigned int s, const StateType& mu, const StateType & T ) const;

  private:

        template <typename StateType>
        const
        ANTIOCH_AUTO(StateType)
         A(const StateType & rho_times_self_diff_over_mu) const
        ANTIOCH_AUTOFUNC(StateType,five_over_two - rho_times_self_diff_over_mu)

        template <typename StateType>
        const
        ANTIOCH_AUTO(StateType)
         B(unsigned int s, const StateType & T, const StateType & rho_times_self_diff_over_mu) const
        ANTIOCH_AUTOFUNC(StateType,_rot(T) + two_over_pi * (five_over_three * this->_thermo.cv_rot_over_R(s) +
                                                    rho_times_self_diff_over_mu ) )


        /*! never ever use it*/
        KineticsTheoryThermalConductivity();

       //small enough
        RotationalRelaxation<CoeffType> _rot;

//constants
        const CoeffType five_over_two;
        const CoeffType five_over_three;
        const CoeffType two_over_pi;
        const CoeffType one;

  };

  template <typename ThermoEvaluator, typename CoeffType>
  inline
  KineticsTheoryThermalConductivity<ThermoEvaluator,CoeffType>::KineticsTheoryThermalConductivity(const ThermoEvaluator & t, CoeffType Z_298K, CoeffType LJ_depth)
    : SpeciesConductivityBase<KineticsTheoryThermalConductivity<ThermoEvaluator,CoeffType> >(),
        _thermo(t),
        _rot(Z_298K,LJ_depth),
        five_over_two(CoeffType(5)/CoeffType(2)),
        five_over_three(CoeffType(5)/CoeffType(3)),
        two_over_pi(2/Constants::pi<CoeffType>()),
        one(1)
  {}

  template <typename ThermoEvaluator, typename CoeffType>
  inline
  KineticsTheoryThermalConductivity<ThermoEvaluator,CoeffType>::KineticsTheoryThermalConductivity( const ThermoEvaluator& t, const std::vector<CoeffType>& coeffs)
    :SpeciesConductivityBase<KineticsTheoryThermalConductivity<ThermoEvaluator,CoeffType> >(),
     _thermo(t),
     _rot(coeffs[0],coeffs[1]),
     five_over_two(CoeffType(5)/CoeffType(2)),
     five_over_three(CoeffType(5)/CoeffType(3)),
     two_over_pi(2/Constants::pi<CoeffType>()),
     one(1)
  {}

  template <typename ThermoEvaluator, typename CoeffType>
  inline
  void KineticsTheoryThermalConductivity<ThermoEvaluator,CoeffType>::reset_coeffs(const std::vector<CoeffType> & coeffs)
  {
      antioch_assert_equal_to(coeffs.size(),2);

      _rot.reset_coeffs(coeffs[0],coeffs[1]);
  }

  template <typename ThermoEvaluator, typename CoeffType>
  inline
  KineticsTheoryThermalConductivity<ThermoEvaluator,CoeffType>::~KineticsTheoryThermalConductivity()
  {
     return;
  }

  template <typename ThermoEvaluator, typename CoeffType>
  template <typename StateType>
  inline
  const ANTIOCH_AUTO(StateType)
        KineticsTheoryThermalConductivity<ThermoEvaluator,CoeffType>::op_with_diff_impl(const unsigned int s, const StateType& mu_s, const StateType & T, const StateType & rho_s, const StateType & Dss) const
  {
      StateType rho_d_m = rho_s * Dss / mu_s; // only once instead of twice
      StateType A_B = two_over_pi * this->A(rho_d_m) / this->B(s, T, rho_d_m);

      return ( mu_s  * this->_thermo.cv_trans(s) * five_over_two * (one - this->_thermo.cv_rot_over_R(s) / this->_thermo.cv_trans_over_R(s) * A_B) +
               rho_s * Dss  * ( this->_thermo.cv_rot(s) * (one + A_B) +
                              this->_thermo.cv_vib(s,T) ) );
  }

  template <typename ThermoEvaluator, typename CoeffType>
  template <typename StateType>
  ANTIOCH_AUTO(StateType)
    KineticsTheoryThermalConductivity<ThermoEvaluator,CoeffType>::op_no_diff_impl( const unsigned int s, const StateType& mu, const StateType & T ) const
  {
    antioch_error();

    /*The following is dummy*/
    return trans(s,mu,T,mu,mu);
  }

} //end namespace Antioch

#endif
