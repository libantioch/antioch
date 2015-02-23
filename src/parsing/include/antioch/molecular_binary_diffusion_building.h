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
#include "antioch/physical_set.h"

// C++
#include <iostream>
#include <vector>


namespace Antioch
{
  // forward declaration
  template <typename CoeffType>
  class TransportSpecies;

  template <typename ThermoEvaluator,typename CoeffType>
  class TransportMixture;

  template<class NumericType, typename Interpolator, typename ThermoEvaluator>
  void build_molecular_binary_diffusion( PhysicalSet<MolecularBinaryDiffusion<NumericType, Interpolator> , TransportMixture<ThermoEvaluator,NumericType> >& D);

// ----------------------------------------- //

  template<typename NumericType, typename Interpolator, typename ThermoEvaluator>
  void build_molecular_binary_diffusion( PhysicalSet<MolecularBinaryDiffusion<NumericType, Interpolator>, TransportMixture<ThermoEvaluator,NumericType> >& D)
  {
     for(unsigned int s = 0; s < D.mixture().n_species(); s++)
     {
        for (unsigned int j = 0; j <= s; j++)
        {
           Initializer<MolecularBinaryDiffusion<NumericType>, bimolecular_diffusion_tag > 
                init(j,*(D.mixture().transport_species()[s]), *(D.mixture().transport_species()[j]));

           D.add_model(D.mixture().species_inverse_name_map().at(s),init);
        }
     }
  }

}

#endif
