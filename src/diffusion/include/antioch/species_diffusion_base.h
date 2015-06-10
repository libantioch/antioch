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

#ifndef ANTIOCH_SPECIES_DIFFUSION_BASE_H
#define ANTIOCH_SPECIES_DIFFUSION_BASE_H

//C++
#include <vector>

// Antioch
#include "antioch/metaprogramming_decl.h" // ANTIOCH_AUTO*

namespace Antioch
{
  //! Base class for species diffusion models
  /*!
   * Species diffusion models are those that directly compute species
   * diffusion coefficients (in constrast to binary diffusion models
   * that are used to compute binary diffusion matrix). We use the
   * curiously recurring template pattern; subclasses need to implement
   * D_impl - computes species diffusion coefficient
   * reset_coeffs_impl - reset model coefficients
   */
  template<typename Subclass, typename CoeffType>
  class SpeciesDiffusionBase
  {
  public:

    SpeciesDiffusionBase(){};

    virtual ~SpeciesDiffusionBase(){};

    //! Reset model coefficients
    void reset_coeffs( std::vector<CoeffType>& coeffs );

    //! Compute species diffusivity
    template<typename StateType>
    ANTIOCH_AUTO(StateType)
    D( const StateType& rho, const StateType& cp, const StateType& k ) const
    ANTIOCH_AUTOFUNC(StateType, static_cast<const Subclass*>(this)->D_impl(rho,cp,k) )

  };

  template<typename Subclass, typename CoeffType>
  void SpeciesDiffusionBase<Subclass,CoeffType>::reset_coeffs( std::vector<CoeffType>& coeffs )
  {
    static_cast<const Subclass*>(this)->reset_coeffs_impl(coeffs);
  }

} // end namespace Antioch

#endif // ANTIOCH_SPECIES_DIFFUSION_BASE_H
