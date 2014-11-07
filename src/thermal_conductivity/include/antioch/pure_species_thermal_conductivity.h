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

#ifndef ANTIOCH_PURE_SPECIES_THERMAL_CONDUCTIVITY_H
#define ANTIOCH_PURE_SPECIES_THERMAL_CONDUCTIVITY_H

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/metaprogramming_decl.h"
#include "antioch/rotational_relaxation.h"
#include "antioch/physics_metaprogramming_decl.h"

//C++
#include <vector>

namespace Antioch{

  template <typename ThermoEvaluator, typename CoeffType>
  class PureSpeciesThermalConductivity
  {
      public:

        PureSpeciesThermalConductivity(const ThermoEvaluator & t, const CoeffType & Z_298K, const CoeffType & LJ_depth);
                                   
        ~PureSpeciesThermalConductivity();

        void reset_coeffs(const std::vector<CoeffType> & coeffs);

        template <typename StateType>
        const ANTIOCH_AUTO(StateType) 
                operator()(unsigned int s, const StateType& mu, const StateType & T, const StateType & rho, const StateType & Dss) const;

        template <typename StateType>
        ANTIOCH_AUTO(StateType)
        trans( const unsigned int s, const StateType& mu, const StateType & T, const StateType & rho, const StateType & Dss) const
        ANTIOCH_AUTOFUNC(StateType, mu * _thermo.cv_trans(s) * five_over_two * (one - two_over_pi * _thermo.cv_rot_over_R(s) / _thermo.cv_trans_over_R(s) * this->A(rho * Dss / mu) / this->B(s, T, rho * Dss / mu)))

        template <typename StateType>
        ANTIOCH_AUTO(StateType)
        rot( const unsigned int s, const StateType& mu, const StateType & T, const StateType & rho, const StateType & Dss) const
        ANTIOCH_AUTOFUNC(StateType, rho * Dss * _thermo.cv_rot(s) * (one + two_over_pi * this->A(rho * Dss / mu) / this->B(s, T, rho * Dss / mu)))

        template <typename StateType>
        ANTIOCH_AUTO(StateType)
        vib( const unsigned int s, const StateType& mu, const StateType & T, const StateType & rho, const StateType & Dss) const
        ANTIOCH_AUTOFUNC(StateType, rho * Dss * _thermo.cv_vib(s,T))

        //! const ref to the rotational relaxation
        const RotationalRelaxation<CoeffType> & rot() const {return _rot;}

        //! const ref to the thermo
        const ThermoEvaluator & thermo() const {return _thermo;}

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
        ANTIOCH_AUTOFUNC(StateType,_rot(T) + two_over_pi * (five_over_three * _thermo.cv_rot_over_R(s) +
                                                    rho_times_self_diff_over_mu ) )


        /*! never ever use it*/
        PureSpeciesThermalConductivity();

        const ThermoEvaluator & _thermo;

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
  PureSpeciesThermalConductivity<ThermoEvaluator,CoeffType>::PureSpeciesThermalConductivity(const ThermoEvaluator & t, const CoeffType & Z_298K, const CoeffType & LJ_depth):
        _thermo(t),
        _rot(Z_298K,LJ_depth),
        five_over_two(5.L/2.L),
        five_over_three(5.L/3.L),
        two_over_pi(2.L/Constants::pi<CoeffType>()),
        one(1.L)
  {
      return;
  }

  template <typename ThermoEvaluator, typename CoeffType>
  inline
  void PureSpeciesThermalConductivity<ThermoEvaluator,CoeffType>::reset_coeffs(const std::vector<CoeffType> & coeffs)
  {
      antioch_assert_equal_to(coeffs.size(),2);

      _rot.reset_coeffs(coeffs[0],coeffs[1]);
  }

  template <typename ThermoEvaluator, typename CoeffType>
  inline
  PureSpeciesThermalConductivity<ThermoEvaluator,CoeffType>::~PureSpeciesThermalConductivity()
  {
     return;
  }

  template <typename ThermoEvaluator, typename CoeffType>
  template <typename StateType>
  inline
  const ANTIOCH_AUTO(StateType) 
        PureSpeciesThermalConductivity<ThermoEvaluator,CoeffType>::operator()(const unsigned int s, const StateType& mu, const StateType & T, const StateType & rho, const StateType & Dss) const
  {
      StateType rho_d_m = rho * Dss / mu; // only once instead of twice
      StateType A_B = two_over_pi * this->A(rho_d_m) / this->B(s, T, rho_d_m);

      return ( mu  * _thermo.cv_trans(s) * five_over_two * (one - _thermo.cv_rot_over_R(s) / _thermo.cv_trans_over_R(s) * A_B) +
               rho * Dss  * ( _thermo.cv_rot(s) * (one + A_B) +
                              _thermo.cv_vib(s,T) ) );
  }

} //end namespace Antioch

#include "antioch/pure_species_the_con_utils_decl.h"

#endif
