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

// Antioch

// Declare metaprogramming overloads before they're used
#include "antioch/eigen_utils_decl.h"
#include "antioch/metaphysicl_utils_decl.h"
#include "antioch/valarray_utils_decl.h"
#include "antioch/vector_utils_decl.h"

#include "antioch/antioch_asserts.h"
#include "antioch/chemical_species.h"
#include "antioch/chemical_mixture.h"
#include "antioch/reaction_set.h"
#include "antioch/read_reaction_set_data_xml.h"
#include "antioch/cea_thermo.h"
#include "antioch/kinetics_evaluator.h"

#include "antioch/eigen_utils.h"
#include "antioch/metaphysicl_utils.h"
#include "antioch/valarray_utils.h"
#include "antioch/vector_utils.h"

static const unsigned int n_species = 5;

template <typename SpeciesVector1, typename SpeciesVector2>
int species_vec_compare (const SpeciesVector1 &a, const SpeciesVector2 &b, const std::string &name)
{
  int found_errors = 0;

  typedef typename Antioch::value_type<SpeciesVector1>::type StateType;
  typedef typename Antioch::value_type<StateType>::type Scalar;

  for (unsigned int s=0; s != n_species; s++)
    {
      using std::abs;
      using std::max;

      const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;

      // Break this expression up to workaround bugs in my Eigen
      // version - RHS
      const StateType rel_error = (a[s] - b[s])/max(a[s],b[s]);
      const StateType abs_rel_error = abs(rel_error);

      if( Antioch::max(abs_rel_error) > tol )
	{
	  found_errors++;
	}
    }

  if (found_errors)
    {
      std::cerr << "Error: Mismatch in vectors " << name << std::endl;
      for( unsigned int s = 0; s < n_species; s++)
	{
	  std::cout << std::scientific << std::setprecision(16)
		    << "a(" << s << ") = " << a[s]
		    << ", b(" << s << ") = " << b[s]
		    << ", a-b(" << s << ") = " << StateType(a[s]-b[s])
		    << std::endl;
	}
    }

  return found_errors;
}

template <typename PairScalars>
int vectester(const std::string& input_name, const PairScalars& example)
{
  typedef typename Antioch::value_type<PairScalars>::type Scalar;

  std::vector<std::string> species_str_list;
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

  PairScalars T = example;
  T[0] = 1500.0;
  T[1] = 1500.0;
  PairScalars P = example;
  P[0] = 1.0e5;
  P[1] = 1.0e5;

  // Mass fractions
  PairScalars massfrac = example;
  massfrac[0] = 0.2;
  massfrac[1] = 0.2;
  const std::vector<PairScalars> Y(n_species,massfrac);

  const PairScalars R_mix = chem_mixture.R(Y);

  const PairScalars rho = P/(R_mix*T);

  std::vector<PairScalars> molar_densities(n_species, example);
  chem_mixture.molar_densities(rho,Y,molar_densities);

  std::vector<PairScalars> h_RT_minus_s_R(n_species, example);
  typedef typename Antioch::CEAThermodynamics<Scalar>::
    template Cache<PairScalars> Cache;
  thermo.h_RT_minus_s_R(Cache(T),h_RT_minus_s_R);

  Antioch::KineticsEvaluator<Scalar> kinetics( reaction_set );
  std::vector<PairScalars> omega_dot(n_species, example);
  
  kinetics.compute_mass_sources( T, rho, R_mix, Y, molar_densities, h_RT_minus_s_R, omega_dot );

  int return_flag = 0;

#ifdef ANTIOCH_HAVE_EIGEN
  {
    typedef Eigen::Array<PairScalars,n_species,1> SpeciesVecEigenType;
    SpeciesVecEigenType eigen_Y = SpeciesVecEigenType::Constant(massfrac);

    const PairScalars eigen_R = chem_mixture.R(eigen_Y);

    return_flag +=
      species_vec_compare(eigen_R,R_mix,"eigen_R");

    SpeciesVecEigenType eigen_molar_densities;

    chem_mixture.molar_densities(rho,eigen_Y,eigen_molar_densities);

    return_flag +=
      species_vec_compare(eigen_molar_densities, molar_densities,
                          "eigen_molar_densities");

    SpeciesVecEigenType eigen_h_RT_minus_s_R;

    thermo.h_RT_minus_s_R(Cache(T),eigen_h_RT_minus_s_R);

    return_flag +=
      species_vec_compare(eigen_h_RT_minus_s_R, h_RT_minus_s_R,
                          "eigen_h_RT_minus_s_R");

    SpeciesVecEigenType eigen_omega_dot;
  
    kinetics.compute_mass_sources( T, rho, eigen_R, eigen_Y, eigen_molar_densities, eigen_h_RT_minus_s_R, eigen_omega_dot );

    return_flag +=
      species_vec_compare(eigen_omega_dot,omega_dot,"eigen_omega_dot");
  }

  {
    typedef Eigen::Matrix<PairScalars,Eigen::Dynamic,1> SpeciesVecEigenType;
    SpeciesVecEigenType eigen_Y =
      SpeciesVecEigenType::Constant(n_species, 1, massfrac);

    const PairScalars eigen_R = chem_mixture.R(eigen_Y);

    return_flag +=
      species_vec_compare(eigen_R,R_mix,"eigen_R");

    SpeciesVecEigenType eigen_molar_densities(n_species,1);

    chem_mixture.molar_densities(rho,eigen_Y,eigen_molar_densities);

    return_flag +=
      species_vec_compare(eigen_molar_densities, molar_densities,
                          "eigen_molar_densities");

    SpeciesVecEigenType eigen_h_RT_minus_s_R(n_species,1);

    thermo.h_RT_minus_s_R(Cache(T),eigen_h_RT_minus_s_R);

    return_flag +=
      species_vec_compare(eigen_h_RT_minus_s_R, h_RT_minus_s_R,
                          "eigen_h_RT_minus_s_R");

    SpeciesVecEigenType eigen_omega_dot(n_species,1);
  
    kinetics.compute_mass_sources( T, rho, eigen_R, eigen_Y, eigen_molar_densities, eigen_h_RT_minus_s_R, eigen_omega_dot );

    return_flag +=
      species_vec_compare(eigen_omega_dot,omega_dot,"eigen_omega_dot");
  }

#endif // ANTIOCH_HAVE_EIGEN

  for( unsigned int s = 0; s < n_species; s++)
    {
      std::cout << std::scientific << std::setprecision(16)
		<< "omega_dot(" << chem_mixture.chemical_species()[s]->species() << ") = "
		<< omega_dot[s] << std::endl;
    }

  std::vector<PairScalars> omega_dot_reg(n_species,example);
  omega_dot_reg[0][0] =  7.9004530802650654e+04;
  omega_dot_reg[1][0] = -3.4113853617637843e+05;
  omega_dot_reg[2][0] = -1.8898881857838202e+05;
  omega_dot_reg[3][0] =  2.1551399274321867e+05;
  omega_dot_reg[4][0] =  2.3560883120889112e+05;
  omega_dot_reg[0][1] = omega_dot_reg[0][0];
  omega_dot_reg[1][1] = omega_dot_reg[1][0];
  omega_dot_reg[2][1] = omega_dot_reg[2][0];
  omega_dot_reg[3][1] = omega_dot_reg[3][0];
  omega_dot_reg[4][1] = omega_dot_reg[4][0];

  return_flag +=
    species_vec_compare(omega_dot, omega_dot_reg, "omega_dot");

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

  returnval +=
    vectester (argv[1], std::valarray<float>(2));
  returnval +=
    vectester (argv[1], std::valarray<double>(2));
// We're not getting the full long double precision yet?
//  returnval = returnval ||
//    vectester (argv[1], std::valarray<long double>(2));
#ifdef ANTIOCH_HAVE_EIGEN
  returnval +=
    vectester (argv[1], Eigen::Array2f());
  returnval +=
    vectester (argv[1], Eigen::Array2d());
// We're not getting the full long double precision yet?
//  returnval = returnval ||
//    vectester (argv[1], Eigen::Array<long double, 2, 1>());
#endif
#ifdef ANTIOCH_HAVE_METAPHYSICL
  returnval +=
    vectester (argv[1], MetaPhysicL::NumberArray<2, float> (0));
  returnval +=
    vectester (argv[1], MetaPhysicL::NumberArray<2, double> (0));
//  returnval = returnval ||
//    vectester (argv[1], MetaPhysicL::NumberArray<2, long double> (0));
#endif

  std::cout << "Found " << returnval << " errors" << std::endl;

  return (returnval != 0) ? 1 : 0;
}
