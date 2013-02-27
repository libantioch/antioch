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

#include "antioch_config.h"

#include <valarray>

#ifdef ANTIOCH_HAVE_EIGEN
#include "Eigen/Dense"
#endif

#ifdef ANTIOCH_HAVE_METAPHYSICL
#include "metaphysicl/numberarray.h"
#endif

// C++
#include <limits>
#include <string>
#include <vector>

#include "antioch/eigen_utils.h"
#include "antioch/metaphysicl_utils.h"
#include "antioch/valarray_utils.h"

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/chemical_species.h"
#include "antioch/chemical_mixture.h"
#include "antioch/reaction_set.h"
#include "antioch/read_reaction_set_data_xml.h"
#include "antioch/cea_thermo.h"
#include "antioch/kinetics.h"

template <typename Scalar, typename PairScalars>
int vectester(const std::string& input_name, const PairScalars& example)
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
  Antioch::ReactionSet<Scalar> reaction_set( chem_mixture );
  Antioch::CEAThermodynamics<Scalar> thermo( chem_mixture );

  Antioch::read_reaction_set_data_xml<Scalar>( input_name, true, reaction_set );

  Antioch::Kinetics<Scalar> kinetics( reaction_set );
  std::vector<PairScalars> omega_dot(n_species, example);

  PairScalars P = example;
  P[0] = 1.0e5;
  P[1] = 1.2e5;

  // Mass fractions
  PairScalars Y_vals = example;
  Y_vals[0] = 0.2;
  Y_vals[1] = 0.25;
  std::vector<PairScalars> Y(n_species,Y_vals);
  Y[5][1] = 0;

  const PairScalars R_mix = chem_mixture.R(Y);

  const unsigned int n_T_samples = 10;
  const Scalar T0 = 500;
  const Scalar T_inc = 500;

  std::vector<PairScalars> molar_densities(n_species,example);
  std::vector<PairScalars> h_RT_minus_s_R(n_species, example);

  int return_flag = 0;

  for( unsigned int i = 0; i < n_T_samples; i++ )
    {
      PairScalars T = example;
      T[0] = T0 + T_inc*static_cast<Scalar>(i);
      T[1] = T[0]+T_inc/2;

      const PairScalars rho = P/(R_mix*T);
      chem_mixture.molar_densities(rho,Y,molar_densities);

      typedef typename Antioch::CEAThermodynamics<Scalar>::
	template Cache<PairScalars> Cache;

      thermo.h_RT_minus_s_R(Cache(T),h_RT_minus_s_R);

      kinetics.compute_mass_sources( T, rho, R_mix, Y, molar_densities, h_RT_minus_s_R, omega_dot );

      // Omega dot had better sum to 0.0
      PairScalars sum = omega_dot[0];
      for( unsigned int s = 1; s < n_species; s++ )
	{
	  sum += omega_dot[s];
	}
      const Scalar sum_tol = std::numeric_limits<Scalar>::epsilon() * 1e6; // 1.0e-10;
      const PairScalars abs_sum = std::abs(sum);
      if( Antioch::max(abs_sum) > sum_tol )
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

  int returnval = 0;

  returnval = returnval ||
    vectester<float, std::valarray<float> >
      (argv[1], std::valarray<float>(2));
  returnval = returnval ||
    vectester<double, std::valarray<double> >
      (argv[1], std::valarray<double>(2));
  returnval = returnval ||
    vectester<long double, std::valarray<long double> >
      (argv[1], std::valarray<long double>(2));
#ifdef ANTIOCH_HAVE_EIGEN
  returnval = returnval ||
    vectester<float, Eigen::Array2f>
      (argv[1], Eigen::Array2f());
  returnval = returnval ||
    vectester<double, Eigen::Array2d>
      (argv[1], Eigen::Array2d());
  returnval = returnval ||
    vectester<long double, Eigen::Array<long double, 2, 1> >
      (argv[1], Eigen::Array<long double, 2, 1>());
#endif
#ifdef ANTIOCH_HAVE_METAPHYSICL
  returnval = returnval ||
    vectester<float, MetaPhysicL::NumberArray<2, float> >
      (argv[1], 0);
  returnval = returnval ||
    vectester<double, MetaPhysicL::NumberArray<2, double> >
      (argv[1], 0);
  returnval = returnval ||
    vectester<long double, MetaPhysicL::NumberArray<2, long double> >
      (argv[1], 0);
#endif

  return returnval;
}
