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

#ifndef ANTIOCH_SUTHERLAND_VISCOSITY_H
#define ANTIOCH_SUTHERLAND_VISCOSITY_H

// C++
#include <cmath>

namespace Antioch
{
  template<class NumericType>
  class SutherlandViscosity
  {
  public:

    SutherlandViscosity( const NumericType mu_ref, const NumericType T_ref );
    ~SutherlandViscosity();

    NumericType operator()( NumericType T ) const;
    
    void reset_coeffs( const NumericType mu_ref, const NumericType T_ref );

    //! Formatted print, by default to \p std::cout
    void print(std::ostream& os = std::cout) const;

    //! Formatted print.
    friend std::ostream& operator<<(std::ostream& os, const SutherlandViscosity& rate)
    {
      rate.print(os);
      return os;
    }

  protected:

    NumericType _mu_ref;
    NumericType _T_ref;
    
  private:
    
    SutherlandViscosity();

  };

  template<class NumericType>
  SutherlandViscosity<NumericType>::SutherlandViscosity( const NumericType mu_ref, const NumericType T_ref )
    : _mu_ref(mu_ref), _T_ref(T_ref)
  {
    return;
  }

  template<class NumericType>
  SutherlandViscosity<NumericType>::~SutherlandViscosity()
  {
    return;
  }

  template<class NumericType>
  void SutherlandViscosity<NumericType>::print(std::ostream& os) const
  {
    os << _mu_ref << "*T^(3/2)/(T + " << _T_ref << ")" << std::endl;

    return;
  }

  /* ------------------------- Inline Functions ------------------------- */
  template<class NumericType>
  inline
  NumericType SutherlandViscosity<NumericType>::operator()( NumericType T ) const
  {
    return _mu_ref*std::pow(T,1.5)/(T+_T_ref);
  }

  template<class NumericType>
  inline
  void SutherlandViscosity<NumericType>::reset_coeffs( const NumericType mu_ref, const NumericType T_ref )
  {
    _mu_ref = mu_ref;
    _T_ref = T_ref;
    return;
  }

} // end namespace Antioch

#endif //ANTIOCH_SUTHERLAND_VISCOSITY_H
