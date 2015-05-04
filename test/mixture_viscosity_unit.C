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
#include "antioch/mixture_viscosity.h"
#include "antioch/blottner_parsing.h"
#include "antioch/sutherland_parsing.h"

#include "antioch/nasa_mixture.h"
#include "antioch/nasa_evaluator.h"
#include "antioch/nasa_mixture_parsing.h"
#include "antioch/thermo_handler.h"
#include "antioch/stat_mech_thermo.h"

template <typename Scalar>
int tester()
{
  std::vector<std::string> species_str_list;
  const unsigned int n_species = 2;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "Air" ); // Yes, I know this doesn't make sense, it's just a test.

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );

  typedef Antioch::StatMechThermodynamics<Scalar> MicroThermo;

  MicroThermo thermo_stat( chem_mixture );

  Antioch::NASAThermoMixture<Scalar,Antioch::NASA9CurveFit<Scalar> > cea_mixture( chem_mixture );
  Antioch::read_nasa_mixture_data( cea_mixture, Antioch::DefaultFilename::thermo_data(),Antioch::ASCII, true );
  Antioch::NASAEvaluator<Scalar,Antioch::NASA9CurveFit<Scalar> > thermo_mix( cea_mixture );

  typedef Antioch::ThermoHandler<Scalar,Antioch::NASAEvaluator<Scalar,Antioch::NASA9CurveFit<Scalar> >, MicroThermo > Thermo;

  Thermo thermo_handler(thermo_mix,thermo_stat);

  Antioch::TransportMixture<Thermo,Scalar> tran_mixture( chem_mixture, thermo_handler );

  Antioch::MixtureViscosity<Antioch::SutherlandViscosity<Scalar>,Thermo,Scalar>
    s_mu_mixture(tran_mixture);

  Antioch::MixtureViscosity<Antioch::BlottnerViscosity<Scalar>,Thermo,Scalar>
    b_mu_mixture(tran_mixture);

  Antioch::read_sutherland_data_ascii<Thermo,Scalar>( s_mu_mixture, Antioch::DefaultFilename::sutherland_data() );
  Antioch::read_blottner_data_ascii<Thermo,Scalar>( b_mu_mixture, Antioch::DefaultFilename::blottner_data() );

  std::cout << s_mu_mixture << std::endl;
  std::cout << b_mu_mixture << std::endl;

  const Scalar T = 1500.1;

  std::cout << "Blottner:" << std::endl;
  for( unsigned int s = 0; s < n_species; s++ )
    {
      std::cout << "mu(" << species_str_list[s] << ") = " << b_mu_mixture(s, T) << std::endl;
    }

  std::cout << "Sutherland:" << std::endl;
  for( unsigned int s = 0; s < n_species; s++ )
    {
      std::cout << "mu(" << species_str_list[s] << ") = " << s_mu_mixture(s, T) << std::endl;
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
