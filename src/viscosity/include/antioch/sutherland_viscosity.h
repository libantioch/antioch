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

#ifndef ANTIOCH_SUTHERLAND_VISCOSITY_H
#define ANTIOCH_SUTHERLAND_VISCOSITY_H

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/cmath_shims.h"

// C++
#include <cmath>
#include <vector>
#include <iostream>

namespace Antioch
{
  template<typename CoeffType=double>
  class SutherlandViscosity
  {
  protected:

    CoeffType _mu_ref;
    CoeffType _T_ref;
    
  public:

    SutherlandViscosity( const CoeffType mu_ref, const CoeffType T_ref );

    SutherlandViscosity( const std::vector<CoeffType>& coeffs );

    ~SutherlandViscosity() {}

    template <typename StateType>
    ANTIOCH_AUTO(StateType)
    operator()( StateType& T ) const
    ANTIOCH_AUTOFUNC(StateType, _mu_ref*ant_pow(T,CoeffType(1.5))/(T+_T_ref))
    
    void reset_coeffs( const CoeffType mu_ref, const CoeffType T_ref );
    void reset_coeffs( const std::vector<CoeffType> coeffs );

    //! Formatted print, by default to \p std::cout
    void print(std::ostream& os = std::cout) const;

    //! Formatted print.
    friend std::ostream& operator<<(std::ostream& os, const SutherlandViscosity& mu)
    {
      mu.print(os);
      return os;
    }

  private:
    
    SutherlandViscosity();

  };

  template<typename CoeffType>
  SutherlandViscosity<CoeffType>::SutherlandViscosity( const CoeffType mu_ref, const CoeffType T_ref )
    : _mu_ref(mu_ref), _T_ref(T_ref)
  {
  }

  template<typename CoeffType>
  SutherlandViscosity<CoeffType>::SutherlandViscosity( const std::vector<CoeffType>& coeffs )
#ifndef NDEBUG
    : _mu_ref(-1), _T_ref(-1)
#else
    : _mu_ref(coeffs[0]), _T_ref(coeffs[1])
#endif
  {
#ifndef NDEBUG
    antioch_assert_equal_to( coeffs.size(), 2 );
    _mu_ref = coeffs[0];
    _T_ref  = coeffs[1];
#endif
  }

  template<typename CoeffType>
  void SutherlandViscosity<CoeffType>::print(std::ostream& os) const
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
  void SutherlandViscosity<CoeffType>::reset_coeffs( const std::vector<CoeffType> coeffs )
  {
    antioch_assert_equal_to(coeffs.size(), 2);
    _mu_ref = coeffs[0];
    _T_ref = coeffs[1];
  }

} // end namespace Antioch

#include "antioch/sutherland_viscosity_utils_decl.h"

#endif //ANTIOCH_SUTHERLAND_VISCOSITY_H
