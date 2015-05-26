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

#ifndef ANTIOCH_DIFFUSION_TRAITS_H
#define ANTIOCH_DIFFUSION_TRAITS_H

#include "antioch/constant_lewis_diffusivity.h"
#include "antioch/molecular_binary_diffusion.h"

namespace Antioch
{
  template<typename DiffModel, typename CoeffType>
  struct DiffusionTraits;

  template<typename CoeffType>
  struct DiffusionTraits<ConstantLewisDiffusivity<CoeffType>,CoeffType>
  {
    static bool const is_species_diffusion = true;
    static bool const is_binary_diffusion = false;
  };

  template<typename CoeffType, typename Interpolator>
  struct DiffusionTraits<MolecularBinaryDiffusion<CoeffType,Interpolator>,CoeffType>
  {
    static bool const is_species_diffusion = false;
    static bool const is_binary_diffusion = true;
  };

} // end namespace Antioch

#endif // ANTIOCH_DIFFUSION_TRAITS_H
