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
#include "antioch/vector_utils_decl.h"
#include "antioch/sutherland_viscosity.h"
#include "antioch/blottner_viscosity.h"
#include "antioch/kinetics_theory_viscosity.h"
#include "antioch/gsl_spliner.h"
#include "antioch/vector_utils.h"

#include "antioch/mixture_viscosity.h"
#include "antioch/blottner_parsing.h"
#include "antioch/sutherland_parsing.h"
#include "antioch/kinetics_theory_viscosity_building.h"

template <typename Scalar>
int test_values( Scalar test_value, Scalar exact_value, Scalar tol )
{
  int return_flag = 0;
  Scalar error = std::abs(test_value - exact_value);

  if( error > tol )
    {
      std::cout << std::setprecision(16) << std::scientific;

      std::cout << "ERROR: Value exceeds tolerance!" << std::endl
      << "test_value  = " << test_value << std::endl
      << "exact_value = " << exact_value << std::endl
      << "error       = " << error << std::endl
      << "tol         = " << tol << std::endl;

      return_flag = 1;
    }

  return return_flag;
}

template <typename Scalar>
int tester()
{
  std::vector<std::string> species_str_list;
  const unsigned int n_species = 2;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "O2" );

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );

  Antioch::TransportMixture<Scalar> tran_mixture( chem_mixture );

  Antioch::SutherlandViscosity<Scalar> s_N2(1.399306e-06, 1.066667e+02);
  Antioch::SutherlandViscosity<Scalar> s_O2(1.693411e-06, 1.270000e+02);

  Antioch::BlottnerViscosity<Scalar> b_N2(2.68142000000e-02,  3.17783800000e-01, -1.13155513000e+01);
  Antioch::BlottnerViscosity<Scalar> b_O2(4.49290000000e-02, -8.26158000000e-02, -9.20194750000e+00);

  Antioch::KineticsTheoryViscosity<Scalar, Antioch::GSLSpliner> k_N2(97.530, 3.621, 0.0, chem_mixture.M(0)/Antioch::Constants::Avogadro<Scalar>());
  Antioch::KineticsTheoryViscosity<Scalar, Antioch::GSLSpliner> k_O2(107.400, 3.458, 0.0, chem_mixture.M(1)/Antioch::Constants::Avogadro<Scalar>() );

  Antioch::MixtureViscosity<Antioch::SutherlandViscosity<Scalar>,Scalar>
    s_mu_mixture(tran_mixture);

  Antioch::MixtureViscosity<Antioch::BlottnerViscosity<Scalar>,Scalar>
    b_mu_mixture(tran_mixture);

  Antioch::MixtureViscosity<Antioch::KineticsTheoryViscosity<Scalar, Antioch::GSLSpliner>,Scalar>
    k_mu_mixture(tran_mixture);

  Antioch::read_sutherland_data_ascii<Scalar>( s_mu_mixture, Antioch::DefaultFilename::sutherland_data() );
  Antioch::read_blottner_data_ascii<Scalar>( b_mu_mixture, Antioch::DefaultFilename::blottner_data() );
  Antioch::build_kinetics_theory_viscosity<Scalar>(k_mu_mixture);

  std::cout << s_mu_mixture << std::endl;
  std::cout << b_mu_mixture << std::endl;
  std::cout << k_mu_mixture << std::endl;

  const Scalar T = 1500.1;

  std::cout << "Sutherland:" << std::endl;
  for( unsigned int s = 0; s < n_species; s++ )
    {
      std::cout << "mu(" << species_str_list[s] << ") = " << s_mu_mixture(s, T) << std::endl;
    }

  std::cout << "Blottner:" << std::endl;
  for( unsigned int s = 0; s < n_species; s++ )
    {
      std::cout << "mu(" << species_str_list[s] << ") = " << b_mu_mixture(s, T) << std::endl;
    }

  std::cout << "Kinetic Theory:" << std::endl;
  for( unsigned int s = 0; s < n_species; s++ )
    {
      std::cout << "mu(" << species_str_list[s] << ") = " << k_mu_mixture(s, T) << std::endl;
    }

  int return_flag = 0;

  Scalar tol = 2.0*std::numeric_limits<Scalar>::epsilon();
  return_flag = test_values( s_mu_mixture(0, T), s_N2(T), tol ) ||
    test_values( s_mu_mixture(1, T), s_O2(T), tol ) ||
    test_values( b_mu_mixture(0, T), b_N2(T), tol ) ||
    test_values( b_mu_mixture(1, T), b_O2(T), tol ) ||
    test_values( k_mu_mixture(0, T), k_N2(T), tol ) ||
    test_values( k_mu_mixture(1, T), k_O2(T), tol );

  return return_flag;
}

int main()
{
  return (tester<double>() ||
          tester<long double>() ||
          tester<float>());
}
