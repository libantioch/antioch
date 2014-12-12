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

#ifndef ANTIOCH_CONSTANT_LEWIS_DIFFUSIVITY_H
#define ANTIOCH_CONSTANT_LEWIS_DIFFUSIVITY_H

//Antioch
#include "antioch/metaprogramming_decl.h" // ANTIOCH_AUTO*

//C++
#include <vector>

namespace Antioch
{

  template<typename CoeffType=double>
  class ConstantLewisDiffusivity
  {
  protected:

    CoeffType _Le;
    
  public:
    
    ConstantLewisDiffusivity( const CoeffType Le ) : _Le(Le) {}
    ConstantLewisDiffusivity( const std::vector<CoeffType> & coefs );

    ~ConstantLewisDiffusivity() {}

    template <typename StateType>
    ANTIOCH_AUTO(StateType)
    operator()( const StateType& rho, const StateType& cp, const StateType& k ) const
    ANTIOCH_AUTOFUNC(StateType, _Le * k / (rho * cp) )

    template<typename StateType>
    ANTIOCH_AUTO(StateType)
    D( const StateType& rho, const StateType& cp, const StateType& k ) const
    ANTIOCH_AUTOFUNC(StateType, _Le*k/(rho*cp))

  };

    template <typename CoeffType>
    ConstantLewisDiffusivity<CoeffType>::ConstantLewisDiffusivity( const std::vector<CoeffType> & coeffs ):
#ifndef NDEBUG
      _Le(-1)
#else
      _Le(coeffs[0])
#endif
    {
#ifndef NDEBUG
      antioch_assert_equal_to(coeffs.size(),1);
      _Le = coeffs[0];
#endif
    }

} // end namespace Antioch

#include "antioch/constant_lewis_diffusivity_utils_decl.h"

#endif // ANTIOCH_CONSTANT_LEWIS_DIFFUSIVITY_H
