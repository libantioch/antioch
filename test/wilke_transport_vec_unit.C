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
#include "antioch/transport_mixture.h"
#include "antioch/default_filename.h"
#include "antioch/stat_mech_thermo.h"
#include "antioch/cea_mixture.h"
#include "antioch/cea_evaluator.h"
#include "antioch/cea_mixture_ascii_parsing.h"
#include "antioch/eucken_thermal_conductivity.h"
#include "antioch/blottner_viscosity.h"
#include "antioch/constant_lewis_diffusivity.h"
#include "antioch/eucken_thermal_conductivity_utils.h"
#include "antioch/blottner_viscosity_utils.h"
#include "antioch/constant_lewis_diffusivity_utils.h"
#include "antioch/physical_set.h"
#include "antioch/physics_metaprogramming.h"
#include "antioch/wilke_mixture.h"
#include "antioch/wilke_evaluator.h"
#include "antioch/blottner_parsing.h"
#include "antioch/eucken_thermal_conductivity_building.h"
#include "antioch/constant_lewis_diffusivity_building.h"

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

  Antioch::StatMechThermodynamics<Scalar> thermo( chem_mixture );

  Antioch::TransportMixture<Antioch::StatMechThermodynamics<Scalar>,Scalar> tran_mixture( chem_mixture, thermo );

  typedef Antioch::TransportMixture<Antioch::StatMechThermodynamics<Scalar>,Scalar> TransportType;
  typedef Antioch::ChemicalMixture<Scalar>                                          ChemicalType;
  typedef Antioch::StatMechThermodynamics<Scalar>                                   ThermoType;

// thermo for cp (diffusion)
  Antioch::CEAThermoMixture<Scalar> cea_mixture( chem_mixture );
  Antioch::read_cea_mixture_data_ascii( cea_mixture, Antioch::DefaultFilename::thermo_data() );
  Antioch::CEAEvaluator<Scalar> thermo_mix( cea_mixture );

  typedef Antioch::CEAEvaluator<Scalar>  ThermoMixType;

  Antioch::PhysicalSet< Antioch::EuckenThermalConductivity<ThermoType>, TransportType > k( tran_mixture );

  Antioch::PhysicalSet<Antioch::BlottnerViscosity<Scalar>, ChemicalType> mu( chem_mixture );

  Antioch::PhysicalSet<Antioch::ConstantLewisDiffusivity<Scalar>, ChemicalType > D( chem_mixture );

  Antioch::read_blottner_data_ascii( mu, Antioch::DefaultFilename::blottner_data() );

  Antioch::build_constant_lewis_diffusivity<Scalar>( D, 1.4);

  typedef Antioch::PhysicalSet< Antioch::EuckenThermalConductivity<ThermoType>, TransportType > TCType;
  typedef Antioch::PhysicalSet<Antioch::BlottnerViscosity<Scalar>, ChemicalType>                VType;
  typedef Antioch::PhysicalSet<Antioch::ConstantLewisDiffusivity<Scalar>, ChemicalType >        DType;
  

//
  Antioch::WilkeMixture<DType,VType,TCType,TransportType,ThermoMixType,Scalar> wilke_mixture(D, mu, k, tran_mixture, thermo_mix );

  typedef Antioch::WilkeMixture<DType,VType,TCType,TransportType,ThermoMixType,Scalar> WilkeType;

  Antioch::WilkeEvaluator< WilkeType > wilke( wilke_mixture);

  int return_flag = 0;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 2;

  // First check phi calculation
  // phi_s = sum_r (chi_r*(1+sqrt(mu_s/mu_r)*(Mr/Ms)^(1/4))^2)/sqrt(8*(1+Ms/Mr))
  {
    std::vector<PairScalars> mu(5, example);
    std::vector<PairScalars> chi(5, example);

    for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
      {
        mu[0][2*tuple  ] = 0.1L;
        mu[1][2*tuple  ] = 0.2L;
        mu[2][2*tuple  ] = 0.3L;
        mu[3][2*tuple  ] = 0.15L;
        mu[4][2*tuple  ] = 0.25L;
        mu[0][2*tuple+1] = 0.25L;
        mu[1][2*tuple+1] = 0.15L;
        mu[2][2*tuple+1] = 0.3L;
        mu[3][2*tuple+1] = 0.2L;
        mu[4][2*tuple+1] = 0.1L;
    
        chi[0][2*tuple  ] = 0.1L;
        chi[1][2*tuple  ] = 0.2L;
        chi[2][2*tuple  ] = 0.3L;
        chi[3][2*tuple  ] = 0.15L;
        chi[4][2*tuple  ] = 0.25L;
        chi[0][2*tuple+1] = 0.25L;
        chi[1][2*tuple+1] = 0.15L;
        chi[2][2*tuple+1] = 0.3L;
        chi[3][2*tuple+1] = 0.2L;
        chi[4][2*tuple+1] = 0.1L;
      }

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

/*   for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
    {
     each_mass[2*tuple  ] = 0.2L;
     each_mass[2*tuple+1] = 0.2L;
    }

   std::vector<PairScalars> mass_fractions( 5, each_mass); 

  // Currently dummy
  const Scalar mu_exact = ;

   PairScalars T = example;
   T[0] = 1000.0L;
   T[1] = 1200.0L;

   PairScalars wilke_mu = wilke.mu(T, mass_fractions );
   PairScalars wilke_k = wilke.k(T, mass_fractions );
*/
  int return_flag_temp = 0;
//  return_flag_temp = test_mu( wilke.mu(T, mass_fractions ), mu_exact, tol );
  if( return_flag_temp != 0 ) return_flag = 1;

  return return_flag;
}

int main()
{
  int returnval = 0;

  returnval = returnval ||
    tester (std::valarray<float>(2*ANTIOCH_N_TUPLES), "valarray<float>");
  returnval = returnval ||
    tester (std::valarray<double>(2*ANTIOCH_N_TUPLES), "valarray<double>");
  returnval = returnval ||
    tester (std::valarray<long double>(2*ANTIOCH_N_TUPLES), "valarray<ld>");
#ifdef ANTIOCH_HAVE_EIGEN
  returnval = returnval ||
    tester (Eigen::Array<float, 2*ANTIOCH_N_TUPLES, 1>(), "Eigen::ArrayXf");
  returnval = returnval ||
    tester (Eigen::Array<double, 2*ANTIOCH_N_TUPLES, 1>(), "Eigen::ArrayXd");
  returnval = returnval ||
    tester (Eigen::Array<long double, 2*ANTIOCH_N_TUPLES, 1>(), "Eigen::ArrayXld");
#endif
#ifdef ANTIOCH_HAVE_METAPHYSICL
  returnval = returnval ||
    tester (MetaPhysicL::NumberArray<2*ANTIOCH_N_TUPLES, float> (0), "NumberArray<float>");
  returnval = returnval ||
    tester (MetaPhysicL::NumberArray<2*ANTIOCH_N_TUPLES, double> (0), "NumberArray<double>");
  returnval = returnval ||
    tester (MetaPhysicL::NumberArray<2*ANTIOCH_N_TUPLES, long double> (0), "NumberArray<ld>");
#endif
#ifdef ANTIOCH_HAVE_VEXCL
  vex::Context ctx_f (vex::Filter::All);
  if (!ctx_f.empty())
    returnval = returnval ||
      tester (vex::vector<float> (ctx_f, 2*ANTIOCH_N_TUPLES), "vex::vector<float>");

  vex::Context ctx_d (vex::Filter::DoublePrecision);
  if (!ctx_d.empty())
    returnval = returnval ||
      tester (vex::vector<double> (ctx_d, 2*ANTIOCH_N_TUPLES), "vex::vector<double>");
#endif

#ifdef ANTIOCH_HAVE_GRVY
  gt.Finalize();
  gt.Summarize();
#endif

  return returnval;
}
