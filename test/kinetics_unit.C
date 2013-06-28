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
#include <limits>
#include <string>
#include <vector>

// Antioch
#include "antioch/vector_utils.h"

#include "antioch/antioch_asserts.h"
#include "antioch/chemical_species.h"
#include "antioch/chemical_mixture.h"
#include "antioch/reaction_set.h"
#include "antioch/read_reaction_set_data_xml.h"
#include "antioch/cea_thermo.h"
#include "antioch/kinetics_evaluator.h"

template <typename Scalar>
int tester_N2N(const std::string& input_name)
{
  using std::abs;

  std::vector<std::string> species_str_list;
  const unsigned int n_species = 2;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "N" );

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );
  Antioch::ReactionSet<Scalar> reaction_set( chem_mixture );
  Antioch::CEAThermodynamics<Scalar> thermo( chem_mixture );

  Antioch::read_reaction_set_data_xml<Scalar>( input_name, false, reaction_set );

  Antioch::KineticsEvaluator<Scalar> kinetics( reaction_set, 0 );
  std::vector<Scalar> omega_dot(n_species);

  const Scalar P = 1.0e5;

  // Mass fractions
  std::vector<Scalar> Y(n_species,0.2);

  const Scalar R_mix = chem_mixture.R(Y);

  const unsigned int n_T_samples = 10;
  const Scalar T0 = 500;
  const Scalar T_inc = 500;

  std::vector<Scalar> molar_densities(n_species,0.0);
  std::vector<Scalar> h_RT_minus_s_R(n_species);

  int return_flag = 0;

  for( unsigned int i = 0; i < n_T_samples; i++ )
    {
      const Scalar T = T0 + T_inc*static_cast<Scalar>(i);
      const Scalar rho = P/(R_mix*T);
      chem_mixture.molar_densities(rho,Y,molar_densities);

      typedef typename Antioch::CEAThermodynamics<Scalar>::template Cache<Scalar> Cache;

      thermo.h_RT_minus_s_R(Cache(T),h_RT_minus_s_R);

      if( i == 0 )
        {
          std::vector<std::vector<Scalar> > loss_matrix;
          std::vector<std::vector<Scalar> > prod_matrix;
          std::vector<std::vector<Scalar> > net_matrix;

          reaction_set.print_chemical_scheme( std::cout, T, molar_densities, h_RT_minus_s_R, loss_matrix, prod_matrix, net_matrix );
        }
      
      kinetics.compute_mass_sources( T, molar_densities, h_RT_minus_s_R, omega_dot );

      // Omega dot had better sum to 0.0
      Scalar sum = 0;
      for( unsigned int s = 0; s < n_species; s++ )
	{
	  sum += omega_dot[s];
	}
      const Scalar sum_tol = std::numeric_limits<Scalar>::epsilon() * 1.0e6; // 1.0e-10;
      if( abs( sum ) > sum_tol )
	{
	  return_flag = 1;
	  std::cerr << "Error: omega_dot did not sum to 0.0." << std::endl
		    << std::scientific << std::setprecision(16)
		    << "T = " << T << std::endl
		    << "sum = " << sum << ", sum_tol = " << sum_tol << std::endl;
	  for( unsigned int s = 0; s < n_species; s++ )
	    {
	      std::cerr << std::scientific << std::setprecision(16)
			<< "omega_dot(" << chem_mixture.chemical_species()[s]->species() << ") = "
			<< omega_dot[s] << std::endl;
	    }
	  std::cout << std::endl << std::endl;
	}
    }
  
  return return_flag;
}


template <typename Scalar>
int tester(const std::string& input_name)
{
  using std::abs;

  std::vector<std::string> species_str_list;
  const unsigned int n_species = 5;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "O2" );
  species_str_list.push_back( "N" );
  species_str_list.push_back( "O" );
  species_str_list.push_back( "NO" );

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );
  Antioch::ReactionSet<Scalar> reaction_set( chem_mixture );
  Antioch::CEAThermodynamics<Scalar> thermo( chem_mixture );

  Antioch::read_reaction_set_data_xml<Scalar>( input_name, false, reaction_set );

  Antioch::KineticsEvaluator<Scalar> kinetics( reaction_set, 0 );
  std::vector<Scalar> omega_dot(n_species);

  const Scalar P = 1.0e5;

  // Mass fractions
  std::vector<Scalar> Y(n_species,0.2);

  const Scalar R_mix = chem_mixture.R(Y);

  const unsigned int n_T_samples = 10;
  const Scalar T0 = 500;
  const Scalar T_inc = 500;

  std::vector<Scalar> molar_densities(n_species,0.0);
  std::vector<Scalar> h_RT_minus_s_R(n_species);

  int return_flag = 0;

  for( unsigned int i = 0; i < n_T_samples; i++ )
    {
      const Scalar T = T0 + T_inc*static_cast<Scalar>(i);
      const Scalar rho = P/(R_mix*T);
      chem_mixture.molar_densities(rho,Y,molar_densities);

      typedef typename Antioch::CEAThermodynamics<Scalar>::template Cache<Scalar> Cache;

      thermo.h_RT_minus_s_R(Cache(T),h_RT_minus_s_R);

      kinetics.compute_mass_sources( T, molar_densities, h_RT_minus_s_R, omega_dot );

      // Omega dot had better sum to 0.0
      Scalar sum = 0;
      for( unsigned int s = 0; s < n_species; s++ )
	{
	  sum += omega_dot[s];
	}
      const Scalar sum_tol = std::numeric_limits<Scalar>::epsilon() * 2.7e6; // 1.6e-10;
      if( abs( sum ) > sum_tol )
	{
	  return_flag = 1;
	  std::cerr << "Error: omega_dot did not sum to 0.0." << std::endl
		    << std::scientific << std::setprecision(16)
		    << "T = " << T << std::endl
		    << "sum = " << sum << ", sum_tol = " << sum_tol << std::endl;
	  for( unsigned int s = 0; s < n_species; s++ )
	    {
	      std::cerr << std::scientific << std::setprecision(16)
			<< "omega_dot(" << chem_mixture.chemical_species()[s]->species() << ") = "
			<< omega_dot[s] << std::endl;
	    }
	  std::cout << std::endl << std::endl;
	}
    }
  
  return return_flag;
}


int main(int argc, char* argv[])
{
  // Check command line count.
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify reaction set XML input file." << std::endl;
      antioch_error();
    }

  int return_flag=0;

  return_flag += (tester<double>(std::string(argv[1])) ||
                  tester<long double>(std::string(argv[1])) ||
                  tester<float>(std::string(argv[1])));

  return_flag += (tester_N2N<double>(std::string(argv[1])) ||
                  tester_N2N<long double>(std::string(argv[1])) ||
                  tester_N2N<float>(std::string(argv[1])));

  return return_flag;
}

