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

#ifndef ANTIOCH_MOLECULAR_BINARY_DIFFUSION_BUILDING_H
#define ANTIOCH_MOLECULAR_BINARY_DIFFUSION_BUILDING_H

// Antioch
#include "antioch/molecular_binary_diffusion.h"
#include "antioch/mixture_binary_diffusion.h"

// C++
#include <iostream>
#include <vector>


namespace Antioch
{
  // forward declaration
  template <typename CoeffType>
  class TransportSpecies;

  template <typename CoeffType>
  class TransportMixture;

  template<typename NumericType, typename Interpolator>
  void build_molecular_binary_diffusion( MixtureBinaryDiffusion<MolecularBinaryDiffusion<NumericType,Interpolator>,NumericType>& D )
  {
    for(unsigned int i = 0; i < D.mixture().n_species(); i++)
     {
        for (unsigned int j = 0; j < D.mixture().n_species(); j++)
        {
           D.add( i, j,
                  D.mixture().transport_species(i),
                  D.mixture().transport_species(j) );
        }
     }
  }

}

#endif
