//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// Antioch - A Gas Dynamics Thermochemistry Library
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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// C++
#include <iostream>

// Antioch
#include "antioch/chemical_mixture.h"
#include "antioch/stat_mech_thermo.h"
#include "antioch/eucken_thermal_conductivity.h"
#include "antioch/blottner_viscosity.h"
#include "antioch/wilke_mixture.h"
#include "antioch/wilke_evaluator.h"
#include "antioch/blottner_parsing.h"

template <typename Scalar>
int test_mu( const Scalar mu, const Scalar mu_exact, const Scalar tol )
{
  int return_flag = 0;

  const Scalar rel_error = std::fabs( (mu - mu_exact)/mu_exact);

  if( rel_error  > tol )
    {
      std::cerr << "Error: Mismatch in thermal conductivity" << std::endl
		<< "mu       = " << mu << std::endl
		<< "mu_exact = " << mu_exact << std::endl
		<< "rel_error = " << rel_error << std::endl
		<< "tol = " << tol << std::endl;
      return_flag = 1;
    }

  return return_flag;
}

template <typename Scalar>
int tester()
{
  std::vector<std::string> species_str_list;
  const unsigned int n_species = 5;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "O2" );
  species_str_list.push_back( "N" );
  species_str_list.push_back( "O" );
  species_str_list.push_back( "NO" );

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );

  const Scalar R_N2 = Antioch::Constants::R_universal<Scalar>()/28.016L;
  const Scalar R_O2 = Antioch::Constants::R_universal<Scalar>()/32.0L;
  const Scalar R_N = Antioch::Constants::R_universal<Scalar>()/14.008L;
  const Scalar R_O = Antioch::Constants::R_universal<Scalar>()/16.0L;
  const Scalar R_NO = Antioch::Constants::R_universal<Scalar>()/30.008L;

  Antioch::WilkeMixture<Scalar> wilke_mixture( chem_mixture );
  
  Antioch::StatMechThermodynamics<Scalar> thermo( chem_mixture );

  Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<Scalar> > k( thermo );

  Antioch::MixtureViscosity<Antioch::BlottnerViscosity<Scalar>, Scalar> mu( chem_mixture );

  Antioch::read_blottner_data_ascii_default( mu );

  Antioch::WilkeEvaluator< Antioch::MixtureViscosity<Antioch::BlottnerViscosity<Scalar>, Scalar>,
                           Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<Scalar> >,
                           Scalar > wilke( wilke_mixture, mu, k );

  const Scalar T = 1000.0L;

  std::vector<Scalar> mass_fractions( 5, 0.2L);

  int return_flag = 0;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 2;

  // Currently dummy
  const Scalar mu_exact = 1.0L;

  wilke.mu(T, mass_fractions );

  int return_flag_temp = 0;
  //return_flag_temp = test_mu( wilke.mu(T, mass_fractions ), mu_exact, tol );
  if( return_flag_temp != 0 ) return_flag = 1;

  return return_flag;
}

int main()
{
  return ( tester<double>() ||
           tester<long double>() ||
           tester<float>() );
}
