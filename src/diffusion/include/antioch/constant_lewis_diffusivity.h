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

#ifndef ANTIOCH_CONSTANT_LEWIS_DIFFUSIVITY_H
#define ANTIOCH_CONSTANT_LEWIS_DIFFUSIVITY_H

// Antioch
#include "antioch/metaprogramming_decl.h" // ANTIOCH_AUTO*
#include "antioch/species_diffusion_base.h"

namespace Antioch
{

  //! Compute species diffusivity based on constant Lewis number
  /*!
   * \f$Le = \frac{k}{\rho c_p D}\f$
   */
  template<typename CoeffType=double>
  class ConstantLewisDiffusivity : public SpeciesDiffusionBase<ConstantLewisDiffusivity<CoeffType>,CoeffType>
  {
  public:

    ConstantLewisDiffusivity( const CoeffType Le );

    ConstantLewisDiffusivity( const std::vector<CoeffType>& coeffs );

    virtual ~ConstantLewisDiffusivity(){}

    //! Friend the base class so we can make the implementation protected
    friend class SpeciesDiffusionBase<ConstantLewisDiffusivity<CoeffType>,CoeffType>;

  protected:
    CoeffType _Le;

    void reset_coeffs_impl( const std::vector<CoeffType>& coeffs );

    template<typename StateType>
    ANTIOCH_AUTO(StateType)
    D_impl( const StateType& rho, const StateType& cp, const StateType& k ) const
    ANTIOCH_AUTOFUNC(StateType, k/(_Le*rho*cp))
  };

  template<typename CoeffType>
  ConstantLewisDiffusivity<CoeffType>::ConstantLewisDiffusivity( const CoeffType Le )
    : SpeciesDiffusionBase<ConstantLewisDiffusivity<CoeffType>,CoeffType>(),
    _Le(Le)
  {}

  template<typename CoeffType>
  ConstantLewisDiffusivity<CoeffType>::ConstantLewisDiffusivity( const std::vector<CoeffType>& coeffs )
    : SpeciesDiffusionBase<ConstantLewisDiffusivity<CoeffType>,CoeffType>(),
    _Le(coeffs[0])
  { antioch_assert_equal_to( coeffs.size(), 1 ); }

  template<typename CoeffType>
  void ConstantLewisDiffusivity<CoeffType>::reset_coeffs_impl( const std::vector<CoeffType>& coeffs )
  {
    antioch_assert_equal_to( coeffs.size(), 1 );
    _Le = coeffs[0];
  }

} // end namespace Antioch

#endif // ANTIOCH_CONSTANT_LEWIS_DIFFUSIVITY_H
