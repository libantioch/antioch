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

#include "antioch/metaprogramming_decl.h" // ANTIOCH_AUTO*

namespace Antioch
{

  template<typename CoeffType=double>
  class ConstantLewisDiffusivity
  {
  protected:

    CoeffType _Le;
    
  public:
    
    ConstantLewisDiffusivity( const CoeffType Le ) : _Le(Le) {}

    ~ConstantLewisDiffusivity() {}
    
    template<typename StateType>
    ANTIOCH_AUTO(StateType)
    D( const StateType& rho, const StateType& cp, const StateType& k ) const
    ANTIOCH_AUTOFUNC(StateType, _Le*k/(rho*cp))

  };

} // end namespace Antioch

#include "antioch/constant_lewis_diffusivity_utils_decl.h"
#include "antioch/constant_lewis_diffusivity_utils.h"

#endif // ANTIOCH_CONSTANT_LEWIS_DIFFUSIVITY_H
