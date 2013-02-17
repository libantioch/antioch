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

#ifndef ANTIOCH_BLOTTNER_VISCOSITY_H
#define ANTIOCH_BLOTTNER_VISCOSITY_H

// C++
#include <cmath>

namespace Antioch
{
  template<class NumericType>
  class BlottnerViscosity
  {
  public:

    BlottnerViscosity( const NumericType a, const NumericType b, const NumericType c );
    ~BlottnerViscosity();

    NumericType operator()( NumericType T ) const;
    
    void reset_coeffs( const NumericType a, const NumericType b, const NumericType c );

    //! Formatted print, by default to \p std::cout
    void print(std::ostream& os = std::cout) const;

    //! Formatted print.
    friend std::ostream& operator<<(std::ostream& os, const BlottnerViscosity& rate)
    {
      rate.print(os);
      return os;
    }

  protected:

    NumericType _a;
    NumericType _b;
    NumericType _c;
    
  private:
    
    BlottnerViscosity();

  };

  template<class NumericType>
  BlottnerViscosity<NumericType>::BlottnerViscosity(const NumericType a,
						    const NumericType b,
						    const NumericType c)
    : _a(a), _b(b), _c(c)
  {
    return;
  }

  template<class NumericType>
  BlottnerViscosity<NumericType>::~BlottnerViscosity()
  {
    return;
  }

  template<class NumericType>
  void BlottnerViscosity<NumericType>::print(std::ostream& os) const
  {
    os << 0.1 << "*exp(" << _a << "*(logT)^2 + " << _b << "*logT + " << _c << ")" << std::endl;

    return;
  }

  /* ------------------------- Inline Functions -------------------------*/
  template<class NumericType>
  inline
  NumericType BlottnerViscosity<NumericType>::operator()( NumericType T ) const
  {
    NumericType logT = std::log(T);
    return 0.1*std::exp( (_a*logT + _b)*logT + _c );
  }

  template<class NumericType>
  inline
  void BlottnerViscosity<NumericType>::reset_coeffs( const NumericType a,
						     const NumericType b,
						     const NumericType c )
  {
    _a = a;
    _b = b;
    _c = c;
    return;
  }

} // end namespace Antioch

#endif //ANTIOCH_BLOTTNER_VISCOSITY_H
