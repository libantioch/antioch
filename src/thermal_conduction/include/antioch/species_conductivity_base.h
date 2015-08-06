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

#ifndef ANTIOCH_SPECIES_CONDUCTIVITY_BASE_H
#define ANTIOCH_SPECIES_CONDUCTIVITY_BASE_H

// Antioch
#include "antioch/metaprogramming_decl.h" // ANTIOCH_AUTO*

namespace Antioch
{
  //! Base class for species conducitivity models
  /*!
   * Defines interface to be used by MixtureConductivity to evaluate
   * species conductivity models. Some models require species diffusion
   * to evaluate the conductivity; these are "with diffusion" models.
   * Other models do not require diffusion; these are "no diffusion" models.
   * We use the curiously recurring template pattern; subclasses should implement
   * op_no_diff_impl - "no diffusion" version (should call antioch_error if invalid)
   * op_with_diff_impl - "requires diffusion" version (should call antioch_error if invalid)
   */
  template<class Subclass>
  class SpeciesConductivityBase
  {
  public:

    SpeciesConductivityBase(){}

    virtual ~SpeciesConductivityBase(){}

    //! Compute species conductivity, without species diffusion
    template <typename StateType>
    ANTIOCH_AUTO(StateType)
      operator()( const unsigned int s, const StateType& mu, const StateType & T ) const
    ANTIOCH_AUTOFUNC(StateType, static_cast<const Subclass*>(this)->op_no_diff_impl(s,mu,T))


    //! Compute species conductivity, with species diffusion
    template <typename StateType>
    const ANTIOCH_AUTO(StateType)
      operator()(unsigned int s, const StateType& mu_s, const StateType & T, const StateType & rho_s, const StateType & Dss) const
    ANTIOCH_AUTOFUNC(StateType, static_cast<const Subclass*>(this)->op_with_diff_impl(s,mu_s,T,rho_s,Dss))

  };
} // end namespace Antioch

#endif // ANTIOCH_SPECIES_CONDUCTIVITY_BASE_H
