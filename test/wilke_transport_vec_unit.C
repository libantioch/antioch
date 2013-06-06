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

#ifdef ANTIOCH_HAVE_EIGEN
#include "Eigen/Dense"
#endif

#ifdef ANTIOCH_HAVE_METAPHYSICL
#include "metaphysicl/numberarray.h"
#endif

#ifdef ANTIOCH_HAVE_VEXCL
#include "vexcl/vexcl.hpp"
#endif

// Antioch
// Declare metaprogramming overloads before they're used
#include "antioch/eigen_utils_decl.h"
#include "antioch/metaphysicl_utils_decl.h"
#include "antioch/valarray_utils_decl.h"
#include "antioch/vector_utils_decl.h"
#include "antioch/vexcl_utils_decl.h"

#include "antioch/chemical_mixture.h"
#include "antioch/stat_mech_thermo.h"
#include "antioch/eucken_thermal_conductivity.h"
#include "antioch/blottner_viscosity.h"
#include "antioch/wilke_mixture.h"
#include "antioch/wilke_evaluator.h"
#include "antioch/blottner_parsing.h"

#include "antioch/eigen_utils.h"
#include "antioch/metaphysicl_utils.h"
#include "antioch/valarray_utils.h"
#include "antioch/vector_utils.h"
#include "antioch/vexcl_utils.h"

#ifdef ANTIOCH_HAVE_GRVY
#include "grvy.h"

GRVY::GRVY_Timer_Class gt;
#endif

// C++
#include <iostream>
#include <cmath>

template <typename Scalar, typename PairScalars>
int test_val( const PairScalars val, const PairScalars val_exact,
              const Scalar tol, const std::string& val_name )
{
  using std::abs;

  int return_flag = 0;

  const PairScalars rel_error = (val - val_exact)/val_exact;
  const PairScalars abs_rel_error = abs(rel_error);

  if( Antioch::max(abs_rel_error)  > tol )
    {
      std::cerr << "Error: Mismatch in " << val_name << std::endl
		<< val_name << "    = " << val << std::endl
		<< val_name+"_exact = " << val_exact << std::endl
		<< "abs_rel_error = " << abs_rel_error << std::endl
		<< "tol = " << tol << std::endl;
      return_flag = 1;
    }

  return return_flag;
}

template <typename PairScalars>
int tester(const PairScalars& example, const std::string& testname)
{
  using std::pow;
  using std::sqrt;

  typedef typename Antioch::value_type<PairScalars>::type Scalar;

  std::vector<std::string> species_str_list;
  const unsigned int n_species = 5;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "O2" );
  species_str_list.push_back( "N" );
  species_str_list.push_back( "O" );
  species_str_list.push_back( "NO" );

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );

  Antioch::WilkeMixture<Scalar> wilke_mixture( chem_mixture );
  
  Antioch::StatMechThermodynamics<Scalar> thermo( chem_mixture );

  Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<Scalar> > k( thermo );

  Antioch::MixtureViscosity<Antioch::BlottnerViscosity<Scalar>, Scalar> mu( chem_mixture );

  Antioch::read_blottner_data_ascii_default( mu );

  Antioch::WilkeEvaluator< Antioch::MixtureViscosity<Antioch::BlottnerViscosity<Scalar>, Scalar>,
                           Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<Scalar> >,
                           Scalar > wilke( wilke_mixture, mu, k );

  int return_flag = 0;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 2;

  // First check phi calculation
  // phi_s = sum_r (chi_r*(1+sqrt(mu_s/mu_r)*(Mr/Ms)^(1/4))^2)/sqrt(8*(1+Ms/Mr))
  {
    std::vector<PairScalars> mu(5, example);
    mu[0][0] = 0.1L;
    mu[1][0] = 0.2L;
    mu[2][0] = 0.3L;
    mu[3][0] = 0.15L;
    mu[4][0] = 0.25L;
    mu[0][1] = 0.25L;
    mu[1][1] = 0.15L;
    mu[2][1] = 0.3L;
    mu[3][1] = 0.2L;
    mu[4][1] = 0.1L;

    std::vector<PairScalars> chi(5, example);
    chi[0][0] = 0.1L;
    chi[1][0] = 0.2L;
    chi[2][0] = 0.3L;
    chi[3][0] = 0.15L;
    chi[4][0] = 0.25L;
    chi[0][1] = 0.25L;
    chi[1][1] = 0.15L;
    chi[2][1] = 0.3L;
    chi[3][1] = 0.2L;
    chi[4][1] = 0.1L;

    PairScalars phi_N_exact = Antioch::zero_clone(example);
    unsigned int N_index = 2;
    const Scalar M_N = chem_mixture.M(N_index);

    for( unsigned int r = 0; r < 5; r++ )
      {
        Scalar M_r = chem_mixture.M(r);
        PairScalars dummy = Scalar(1) + sqrt(mu[N_index]/mu[r])*pow( M_r/M_N, Scalar(0.25L) );
        phi_N_exact += chi[r]*dummy*dummy/sqrt(Scalar(8)*( Scalar(1) + M_N/M_r ) );
      }

#ifdef ANTIOCH_HAVE_GRVY
  gt.BeginTimer(testname);
#endif

    const PairScalars phi_N = wilke.compute_phi( mu, chi, N_index );

#ifdef ANTIOCH_HAVE_GRVY
  gt.EndTimer(testname);
#endif

    std::cout << "mu =    " << mu << std::endl;
    std::cout << "chi =   " << chi << std::endl;
    std::cout << "phi_N = " << phi_N << std::endl;

    return_flag = test_val( phi_N, phi_N_exact, tol, std::string("phi") );
  }
  

  PairScalars each_mass = example;
  each_mass[0] = 0.2L;
  each_mass[1] = 0.2L;
  std::vector<PairScalars> mass_fractions( 5, each_mass); 

  // Currently dummy
  //const Scalar mu_exact = ;

  // PairScalars T = example;
  // T[0] = 1000.0L;
  // T[1] = 1200.0L;

  // PairScalars wilke_mu = wilke.mu(T, mass_fractions );
  // PairScalars wilke_k = wilke.k(T, mass_fractions );

  int return_flag_temp = 0;
  //return_flag_temp = test_mu( wilke.mu(T, mass_fractions ), mu_exact, tol );
  if( return_flag_temp != 0 ) return_flag = 1;

  return return_flag;
}

int main()
{
  int returnval = 0;

  returnval = returnval ||
    tester (std::valarray<float>(2), "valarray<float>");
  returnval = returnval ||
    tester (std::valarray<double>(2), "valarray<double>");
  returnval = returnval ||
    tester (std::valarray<long double>(2), "valarray<ld>");
#ifdef ANTIOCH_HAVE_EIGEN
  returnval = returnval ||
    tester (Eigen::Array2f(), "Eigen::Array2f");
  returnval = returnval ||
    tester (Eigen::Array2d(), "Eigen::Array2d");
  returnval = returnval ||
    tester (Eigen::Array<long double, 2, 1>(), "Eigen::Array<ld>");
#endif
#ifdef ANTIOCH_HAVE_METAPHYSICL
  returnval = returnval ||
    tester (MetaPhysicL::NumberArray<2, float> (0), "NumberArray<float>");
  returnval = returnval ||
    tester (MetaPhysicL::NumberArray<2, double> (0), "NumberArray<double>");
  returnval = returnval ||
    tester (MetaPhysicL::NumberArray<2, long double> (0), "NumberArray<ld>");
#endif
#ifdef ANTIOCH_HAVE_VEXCL
  vex::Context ctx (vex::Filter::DoublePrecision);

  returnval = returnval ||
    tester (vex::vector<float> (ctx, 2), "vex::vector<float>");
  returnval = returnval ||
    tester (vex::vector<double> (ctx, 2), "vex::vector<double>");
#endif
  return returnval;
}
