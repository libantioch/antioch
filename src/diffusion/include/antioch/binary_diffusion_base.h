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

#ifndef ANTIOCH_BINARY_DIFFUSION_BASE_H
#define ANTIOCH_BINARY_DIFFUSION_BASE_H

//C++
#include <vector>

// Antioch
#include "antioch/metaprogramming_decl.h" // ANTIOCH_AUTO*

namespace Antioch
{
  template<typename Subclass, typename CoeffType>
  class BinaryDiffusionBase
  {
  public:

    BinaryDiffusionBase(){};

    virtual ~BinaryDiffusionBase(){};

    void reset_coeffs( const TransportSpecies<CoeffType> & s_i,
                       const TransportSpecies<CoeffType> & s_j );

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

} // end namespace Antioch

#endif // ANTIOCH_BINARY_DIFFUSION_BASE_H
