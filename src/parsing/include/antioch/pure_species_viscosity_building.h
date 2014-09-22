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

#ifndef ANTIOCH_PURE_SPECIES_BUILDING_H
#define ANTIOCH_PURE_SPECIES_BUILDING_H

// Antioch
#include "antioch/pure_species_viscosity.h"
#include "antioch/physical_set.h"

// C++
#include <iostream>
#include <vector>


namespace Antioch
{

  template <typename Thermo, typename Scalar>
  class TransportMixture;

  template<class NumericType, class ThermoEvaluator>
  void build_pure_species_viscosity( PhysicalSet<PureSpeciesViscosity<NumericType>, TransportMixture<ThermoEvaluator,NumericType> >& mu);

// ----------------------------------------- //

  template<class NumericType, class ThermoEvaluator>
  void build_pure_species_viscosity( PhysicalSet<PureSpeciesViscosity<NumericType>, TransportMixture<ThermoEvaluator,NumericType> >& mu)
  {
      for(unsigned int s = 0; s < mu.mixture().n_species(); s++)
      {
          std::vector<NumericType> coeffs(4,0);
          coeffs[0] = mu.mixture().transport_species()[s]->LJ_depth();
          coeffs[1] = mu.mixture().transport_species()[s]->LJ_diameter();
          coeffs[2] = mu.mixture().transport_species()[s]->dipole_moment();
          coeffs[3] = mu.mixture().transport_species()[s]->M() / Constants::Avogadro<NumericType>();
          mu.add_model(mu.mixture().species_inverse_name_map().at(s),coeffs);
      }
  }

}

#endif
