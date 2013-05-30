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

// Antioch
#include "antioch/vector_utils.h"

#include "antioch/antioch_asserts.h"
#include "antioch/chemical_species.h"
#include "antioch/chemical_mixture.h"
#include "antioch/reaction_set.h"
#include "antioch/read_reaction_set_data_xml.h"
#include "antioch/cea_thermo.h"
#include "antioch/kinetics_evaluator.h"

// Thrust
#include "thrust/for_each.h"
#include "thrust/device_vector.h"
#include "thrust/iterator/zip_iterator.h"

// C++
#include <limits>
#include <iostream>
#include <string>
#include <vector>

template <typename Scalar>
struct antioch_functor
{
  antioch_functor(const std::vector<std::string> &species_str_list,
		  const std::string& reaction_data_filename) :
    chem_mixture(species_str_list),
    reaction_set(chem_mixture),
    thermo(chem_mixture),
    kinetics(reaction_set, 0)
  { 
    Antioch::read_reaction_set_data_xml<Scalar>
      (reaction_data_filename, false, reaction_set);

    P = 1.0e5;
    Y.resize(species_str_list.size(),1.0/species_str_list.size());
    R_mix = chem_mixture.R(Y);
  }

  template <typename Tuple>
  __host__ __device__
  void operator()(Tuple t)
  {
    const Scalar& T = thrust::get<0>(t);
    const Scalar rho = P / (R_mix * T);

    const short n_species = Y.size();
    std::vector<Scalar> molar_densities(n_species,0.0),
		        h_RT_minus_s_R (n_species),
			omega_dot(n_species);

    chem_mixture.molar_densities(rho,Y,molar_densities);

    typedef typename Antioch::CEAThermodynamics<Scalar>::template Cache<Scalar> Cache;

    thermo.h_RT_minus_s_R(Cache(T),h_RT_minus_s_R);

    kinetics.compute_mass_sources( T, rho, R_mix, Y, molar_densities, h_RT_minus_s_R, omega_dot );
  }

  // Make these first three shared?
  Antioch::ChemicalMixture<Scalar>   chem_mixture;
  Antioch::ReactionSet<Scalar>       reaction_set;
  Antioch::CEAThermodynamics<Scalar> thermo;

  // Need one of these per object
  Antioch::KineticsEvaluator<Scalar> kinetics;

  // Matching the other unit tests
  Scalar P;
  std::vector<Scalar> Y;
  Scalar R_mix;
};

template <typename Scalar>
int tester_N2N(const std::string& input_name)
{
  using std::abs;

  std::vector<std::string> species_str_list;
  const unsigned int n_species = 2;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "N" );

  int return_flag = 0;

/*
  for( unsigned int i = 0; i < n_T_samples; i++ )
    {
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
*/
  
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

  int return_flag = 0;

/*
  for( unsigned int i = 0; i < n_T_samples; i++ )
    {
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
*/
  
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

