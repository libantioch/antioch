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

#ifndef ANTIOCH_WILKE_MIXTURE_H
#define ANTIOCH_WILKE_MIXTURE_H

// Antioch
#include "antioch/wilke_transport_mixture.h"
#include "antioch/physics_placeholder.h"
#include "antioch/antioch_asserts.h"

// C++
#include <vector>

namespace Antioch
{

  // forward declaration to keep the call to 
  // the ChemicalMixture object
  template <typename CoeffType>
  class ChemicalMixture;

  template<class CoeffType = double>  // type
  class WilkeMixture:public WilkeTransportMixture<ChemicalMixture<CoeffType>,PhysicsPlaceholder,CoeffType>
  {
  public:

    WilkeMixture( const ChemicalMixture<CoeffType> & mixture);
    ~WilkeMixture();

  };

  template<class CoeffType>
  WilkeMixture<CoeffType>::WilkeMixture( const ChemicalMixture<CoeffType>& mixture)
    : WilkeTransportMixture<ChemicalMixture<CoeffType>,PhysicsPlaceholder,CoeffType>(mixture,PhysicsPlaceholder())
  {
    antioch_deprecated();
    return;
  }

  template<class CoeffType>
  WilkeMixture<CoeffType>::~WilkeMixture()
  {
    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_WILKE_MIXTURE_H
