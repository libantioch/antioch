//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014-2016 Paul T. Bauman, Benjamin S. Kirk,
//                         Sylvain Plessis, Roy H. Stonger
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

#ifndef ANTIOCH_SUTHERLAND_VISCOSITY_H
#define ANTIOCH_SUTHERLAND_VISCOSITY_H

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/cmath_shims.h"
#include "antioch/species_viscosity_base.h"

// C++
#include <cmath>
#include <vector>
#include <iostream>

namespace Antioch
{
  template<typename CoeffType=double>
  class SutherlandViscosity : public SpeciesViscosityBase<SutherlandViscosity<CoeffType>,CoeffType>
  {
  protected:

    CoeffType _mu_ref;
    CoeffType _T_ref;

  public:

    SutherlandViscosity( const CoeffType mu_ref, const CoeffType T_ref );

    SutherlandViscosity( const std::vector<CoeffType>& coeffs );

    ~SutherlandViscosity(){};

    void reset_coeffs( const CoeffType mu_ref, const CoeffType T_ref );

    //! Friend base class so we can make implementation protected
    friend class SpeciesViscosityBase<SutherlandViscosity<CoeffType>,CoeffType>;

  protected:

    template <typename StateType>
    ANTIOCH_AUTO(StateType)
    op_impl( StateType& T ) const
    ANTIOCH_AUTOFUNC(StateType, _mu_ref*ant_pow(T,CoeffType(1.5))/(T+_T_ref))

    void reset_coeffs_impl( const std::vector<CoeffType>& coeffs );

    void print_impl(std::ostream& os) const;

    //! No extrapolation needed
    /*!
     * Implementation needed for the interface, but we just throw an error
     * is this is called since it's not defined for SutherlandViscosity.
     */
    template <typename StateType>
    void extrapolate_max_temp_impl(const StateType& Tmax);

  private:

    SutherlandViscosity();

  };

  template<typename CoeffType>
  SutherlandViscosity<CoeffType>::SutherlandViscosity( const CoeffType mu_ref, const CoeffType T_ref )
    : SpeciesViscosityBase<SutherlandViscosity<CoeffType>,CoeffType>(),
    _mu_ref(mu_ref), _T_ref(T_ref)
  {
  }

  template<typename CoeffType>
  SutherlandViscosity<CoeffType>::SutherlandViscosity( const std::vector<CoeffType>& coeffs )
#ifndef NDEBUG
    : SpeciesViscosityBase<SutherlandViscosity<CoeffType>,CoeffType>(),
    _mu_ref(-1), _T_ref(-1)
#else
    : SpeciesViscosityBase<SutherlandViscosity<CoeffType>,CoeffType>(),
    _mu_ref(coeffs[0]), _T_ref(coeffs[1])
#endif
  {
#ifndef NDEBUG
    antioch_assert_equal_to( coeffs.size(), 2 );
    _mu_ref = coeffs[0];
    _T_ref  = coeffs[1];
#endif
  }

  template<typename CoeffType>
  void SutherlandViscosity<CoeffType>::print_impl(std::ostream& os) const
  {
    os << _mu_ref << "*T^(3/2)/(T + " << _T_ref << ")" << std::endl;
  }

  template<typename CoeffType>
  inline
  void SutherlandViscosity<CoeffType>::reset_coeffs( const CoeffType mu_ref, const CoeffType T_ref )
  {
    _mu_ref = mu_ref;
    _T_ref = T_ref;
  }

  template<typename CoeffType>
  inline
  void SutherlandViscosity<CoeffType>::reset_coeffs_impl( const std::vector<CoeffType>& coeffs )
  {
    antioch_assert_equal_to(coeffs.size(), 2);
    this->reset_coeffs(coeffs[0], coeffs[1]);
  }

  template<typename CoeffType>
  template <typename StateType>
  inline
  void SutherlandViscosity<CoeffType>::extrapolate_max_temp_impl(const StateType & /*Tmax*/)
  {
    antioch_error_msg("Extrapolation not well defined for SutherlandViscosity!");
  }

} // end namespace Antioch

#endif //ANTIOCH_SUTHERLAND_VISCOSITY_H
