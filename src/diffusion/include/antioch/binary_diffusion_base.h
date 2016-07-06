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

#ifndef ANTIOCH_BINARY_DIFFUSION_BASE_H
#define ANTIOCH_BINARY_DIFFUSION_BASE_H

//C++
#include <vector>

// Antioch
#include "antioch/metaprogramming_decl.h" // ANTIOCH_AUTO*

namespace Antioch
{
  //! Base class for binary diffusion models
  /*!
   * Binary diffusion models are those that are used to compute
   * a binary diffusion matrix (in constrast to species diffusion models
   * that are used to directly compute species diffusivities). We use the
   * curiously recurring template pattern; subclasses need to implement
   * op_impl - computes binary diffusion coefficient for the species s_i, s_j
   * reset_coeffs_impl - reset TransportSpecies s_i, s_j
   */
  template<typename Subclass, typename CoeffType>
  class BinaryDiffusionBase
  {
  public:

    BinaryDiffusionBase(){};

    virtual ~BinaryDiffusionBase(){};

    void reset_coeffs( const TransportSpecies<CoeffType> & s_i,
                       const TransportSpecies<CoeffType> & s_j );

    //! Extrapolate to input maximum temperature, given in [K]
    /*!
     * Some species binary diffusion models, e.g. MolecularBinaryDiffusion, use interpolated
     * quantities for a given temperature range.
     * If the diffusivity is to be evaluated outside that range, an error will occur.
     * This method will reconstruct the interpolation table, but use a linear extrapolation
     * from the max in the existing table to the input maximum temperature.
     *
     * This method is only applicable to a subset of binary diffusion models. Others will
     * throw a runtime error.
     */
    template<typename StateType>
    void extrapolate_max_temp(const StateType& Tmax);

    template <typename StateType>
    const
    ANTIOCH_AUTO(StateType)
    operator()(const StateType & T, const StateType & molar_density) const
    ANTIOCH_AUTOFUNC(StateType, static_cast<const Subclass*>(this)->op_impl(T, molar_density))

  };

  template<typename Subclass, typename CoeffType>
  void BinaryDiffusionBase<Subclass,CoeffType>::reset_coeffs(const TransportSpecies<CoeffType> & s_i,
                                                             const TransportSpecies<CoeffType> & s_j )
  {
    static_cast<const Subclass*>(this)->reset_coeffs_impl(s_i,s_j);
  }

  template<typename Subclass, typename CoeffType>
  template<typename StateType>
  void BinaryDiffusionBase<Subclass,CoeffType>::extrapolate_max_temp(const StateType& Tmax)
  {
    static_cast<Subclass*>(this)->extrapolate_max_temp_impl(Tmax);
  }

} // end namespace Antioch

#endif // ANTIOCH_BINARY_DIFFUSION_BASE_H
