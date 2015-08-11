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

#ifndef ANTIOCH_BLOTTNER_VISCOSITY_H
#define ANTIOCH_BLOTTNER_VISCOSITY_H

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/species_viscosity_base.h"

// C++
#include <cmath>
#include <vector>
#include <iostream>

namespace Antioch
{
  template<typename CoeffType=double>
  class BlottnerViscosity : public SpeciesViscosityBase<BlottnerViscosity<CoeffType>,CoeffType>
  {
  public:

    BlottnerViscosity( const CoeffType a, const CoeffType b, const CoeffType c );

    BlottnerViscosity( const std::vector<CoeffType>& coeffs );

    virtual ~BlottnerViscosity(){};

    void reset_coeffs( const CoeffType a, const CoeffType b, const CoeffType c );

    //! Friend base class so we can make implementation protected
    friend class SpeciesViscosityBase<BlottnerViscosity<CoeffType>,CoeffType>;

  protected:

    template <typename StateType>
    StateType op_impl( const StateType& T ) const;

    void reset_coeffs_impl( const std::vector<CoeffType> coeffs );

    void print_impl(std::ostream& os) const;

    //! No extrapolation needed
    /*!
     * Implementation needed for the interface, but we just throw an error
     * is this is called since it's not defined for BlottnerViscosity.
     */
    template <typename StateType>
    void extrapolate_max_temp_impl(const StateType& Tmax);

    CoeffType _a;
    CoeffType _b;
    CoeffType _c;

  private:

    BlottnerViscosity();

  };

  template<typename CoeffType>
  BlottnerViscosity<CoeffType>::BlottnerViscosity(const CoeffType a,
						  const CoeffType b,
						  const CoeffType c)
    : SpeciesViscosityBase<BlottnerViscosity<CoeffType>,CoeffType>(),
    _a(a), _b(b), _c(c)
  {
    return;
  }

  template<typename CoeffType>
  BlottnerViscosity<CoeffType>::BlottnerViscosity( const std::vector<CoeffType>& coeffs )
    : SpeciesViscosityBase<BlottnerViscosity<CoeffType>,CoeffType>(),
    _a(-1.0), _b(-1.0), _c(-1.0)
  {
    antioch_assert_equal_to( coeffs.size(), 3);
    _a = coeffs[0];
    _b = coeffs[1];
    _c = coeffs[2];
    return;
  }

  template<typename CoeffType>
  void BlottnerViscosity<CoeffType>::print_impl(std::ostream& os) const
  {
    os << 0.1 << "*exp(" << _a << "*(logT)^2 + " << _b << "*logT + " << _c << ")" << std::endl;

    return;
  }

  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType BlottnerViscosity<CoeffType>::op_impl( const StateType& T ) const
  {
    using std::log;
    using std::exp;
    StateType logT = log(T);
    const CoeffType zero_point_one = 0.1L;

    return zero_point_one*exp( (_a*logT + _b)*logT + _c );
  }

  template<typename CoeffType>
  inline
  void BlottnerViscosity<CoeffType>::reset_coeffs( const CoeffType a,
                                                   const CoeffType b,
                                                   const CoeffType c )
  {
    _a = a;
    _b = b;
    _c = c;
    return;
  }

  template<typename CoeffType>
  inline
  void BlottnerViscosity<CoeffType>::reset_coeffs_impl( const std::vector<CoeffType> coeffs )
  {
    antioch_assert_equal_to(coeffs.size(), 3);
    this->reset_coeffs( coeffs[0], coeffs[1], coeffs[2] );
  }

  template<typename CoeffType>
  template <typename StateType>
  inline
  void BlottnerViscosity<CoeffType>::extrapolate_max_temp_impl(const StateType & /*Tmax*/)
  {
    antioch_msg_error("Extrapolation not well defined for BlottnerViscosity!");
    antioch_error();
  }

} // end namespace Antioch

#endif //ANTIOCH_BLOTTNER_VISCOSITY_H
