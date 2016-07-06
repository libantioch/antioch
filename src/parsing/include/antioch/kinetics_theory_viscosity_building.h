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

#include "antioch_config.h"
#ifdef ANTIOCH_HAVE_GSL // if we do not have it, we don't even define the stuff

#ifndef ANTIOCH_KINETICS_THEORY_BUILDING_H
#define ANTIOCH_KINETICS_THEORY_BUILDING_H

// Antioch
#include "antioch/kinetics_theory_viscosity.h"
#include "antioch/mixture_viscosity.h"

// C++
#include <iostream>
#include <vector>


namespace Antioch
{

  template <typename Scalar>
  class TransportMixture;

  template<class NumericType, class SplineType>
  void build_kinetics_theory_viscosity( MixtureViscosity<KineticsTheoryViscosity<NumericType,SplineType>,NumericType >& mu);

  // ----------------------------------------- //

  template<class NumericType, class SplineType>
  void build_kinetics_theory_viscosity( MixtureViscosity<KineticsTheoryViscosity<NumericType,SplineType>,NumericType >& mu)
  {
      for(unsigned int s = 0; s < mu.mixture().n_species(); s++)
      {
          std::vector<NumericType> coeffs(4,0);
          coeffs[0] = mu.mixture().transport_species()[s]->LJ_depth();
          coeffs[1] = mu.mixture().transport_species()[s]->LJ_diameter();
          coeffs[2] = mu.mixture().transport_species()[s]->dipole_moment();
          coeffs[3] = mu.mixture().transport_species()[s]->M() / Constants::Avogadro<NumericType>();
          mu.add(mu.mixture().species_inverse_name_map().at(s),coeffs);
      }
  }

}

#endif // ANTIOCH_KINETICS_THEORY_BUILDING_H

#endif // ANTIOCH_HAVE_GSL
