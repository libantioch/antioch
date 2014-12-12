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

// C++
#include <iostream>
#include <cmath>

// Antioch
#include "antioch_config.h"
#include "antioch/default_filename.h"
#include "antioch/sutherland_viscosity.h"
#include "antioch/blottner_viscosity.h"
#include "antioch/blottner_viscosity_utils.h"
#include "antioch/sutherland_viscosity_utils.h"
#include "antioch/physics_metaprogramming.h"
#include "antioch/physical_set.h"
#include "antioch/blottner_parsing.h"
#include "antioch/sutherland_parsing.h"

template <typename Scalar>
int tester()
{
  std::vector<std::string> species_str_list;
  const unsigned int n_species = 2;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "Air" ); // Yes, I know this doesn't make sense, it's just a test.

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );

  Antioch::PhysicalSet<Antioch::SutherlandViscosity<Scalar>, Antioch::ChemicalMixture<Scalar> > s_mu_mixture(chem_mixture);

  Antioch::PhysicalSet<Antioch::BlottnerViscosity<Scalar>, Antioch::ChemicalMixture<Scalar> > b_mu_mixture(chem_mixture);

  Antioch::read_sutherland_data_ascii<Scalar>( s_mu_mixture, Antioch::DefaultFilename::sutherland_data() );
  Antioch::read_blottner_data_ascii<Scalar>( b_mu_mixture, Antioch::DefaultFilename::blottner_data() );

  std::cout << s_mu_mixture << std::endl;
  std::cout << b_mu_mixture << std::endl;

  const Scalar T = 1500.1;
  Scalar mu(-1);

  std::cout << "Blottner:" << std::endl;
  for( unsigned int s = 0; s < n_species; s++ )
    {
      b_mu_mixture(s, T, mu);
      std::cout << "mu(" << species_str_list[s] << ") = " << mu << std::endl;
    }

  std::cout << "Sutherland:" << std::endl;
  for( unsigned int s = 0; s < n_species; s++ )
    {
      s_mu_mixture(s, T, mu);
      std::cout << "mu(" << species_str_list[s] << ") = " << mu << std::endl;
    }

  int return_flag = 0;

  return return_flag;
}

int main()
{
  return (tester<double>() ||
          tester<long double>() ||
          tester<float>());
}
