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

#ifndef ANTIOCH_CONSTANT_LEWIS_DIFFUSIVITY_BUILDING_H
#define ANTIOCH_CONSTANT_LEWIS_DIFFUSIVITY_BUILDING_H

// Antioch
#include "antioch/constant_lewis_diffusivity.h"
#include "antioch/mixture_diffusion.h"

// C++
#include <iostream>
#include <vector>

namespace Antioch
{

  //
  template <typename NumericType>
  class ChemicalMixture;

  template<class NumericType>
  void build_constant_lewis_diffusivity( MixtureDiffusion<ConstantLewisDiffusivity<NumericType>,NumericType>& D, const NumericType Le )
  {
    std::vector<NumericType> coeffs(1, Le);
    for(unsigned int s = 0; s < D.mixture().n_species(); s++)
      {
        D.add_species_diffusion(s,coeffs);
      }
  }

  template<class NumericType>
  void build_constant_lewis_diffusivity( MixtureDiffusion<ConstantLewisDiffusivity<NumericType>,NumericType>& D, const std::vector<NumericType>& Le )
  {
    for(unsigned int s = 0; s < D.mixture().n_species(); s++)
      {
        std::vector<NumericType> coeffs(1, Le[s]);
        D.add_species_diffusion(s,coeffs);
      }
  }
}

#endif
