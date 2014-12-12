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
#include "antioch/vector_utils_decl.h"

#include "antioch/default_filename.h"
#include "antioch/transport_mixture.h"
#include "antioch/stat_mech_thermo.h"
#include "antioch/cea_mixture.h"
#include "antioch/cea_evaluator.h"
#include "antioch/cea_mixture_ascii_parsing.h"
//
#include "antioch/eucken_thermal_conductivity.h"
#include "antioch/blottner_viscosity.h"
#include "antioch/constant_lewis_diffusivity.h"
#include "antioch/constant_lewis_diffusivity_utils.h"
#include "antioch/blottner_viscosity_utils.h"
#include "antioch/eucken_thermal_conductivity_utils.h"
#include "antioch/pure_species_thermal_conductivity.h"

#ifdef ANTIOCH_HAVE_GSL

#include "antioch/pure_species_viscosity.h"
#include "antioch/molecular_binary_diffusion.h"
#include "antioch/pure_species_viscosity_utils.h"
#include "antioch/molecular_binary_diffusion_utils.h"
#include "antioch/gsl_spliner.h"

#endif

#include "antioch/pure_species_thermal_conductivity_utils.h"
//
#include "antioch/physical_set.h"
#include "antioch/wilke_mixture.h"
#include "antioch/wilke_evaluator.h"
#include "antioch/blottner_parsing.h"
#include "antioch/eucken_thermal_conductivity_building.h"
#include "antioch/constant_lewis_diffusivity_building.h"
#include "antioch/physics_metaprogramming.h"

#include "antioch/vector_utils.h"

template <typename Scalar>
int test_val( const Scalar val, const Scalar val_exact, const Scalar tol, const std::string& val_name )
{
  using std::abs;

  int return_flag = 0;

  const Scalar rel_error = abs( (val - val_exact)/val_exact);

  if( rel_error  > tol )
    {
      std::cerr << "Error: Mismatch in " << val_name << std::endl
		<< val_name << "    = " << val << std::endl
		<< val_name+"_exact = " << val_exact << std::endl
		<< "rel_error = " << rel_error << std::endl
		<< "tol = " << tol << std::endl;
      return_flag = 1;
    }

  return return_flag;
}

template <typename Scalar>
int tester()
{
  using std::pow;

  std::vector<std::string> species_str_list;
  const unsigned int n_species = 5;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "O2" );
  species_str_list.push_back( "N" );
  species_str_list.push_back( "O" );
  species_str_list.push_back( "NO" );

// mixture and thermo for conduction
  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );
  Antioch::StatMechThermodynamics<Scalar> thermo_stat( chem_mixture );
  Antioch::TransportMixture<Antioch::StatMechThermodynamics<Scalar>,Scalar> tran_mixture( chem_mixture, thermo_stat );

// thermo for cp (diffusion)
  Antioch::CEAThermoMixture<Scalar> cea_mixture( chem_mixture );
  Antioch::read_cea_mixture_data_ascii( cea_mixture, Antioch::DefaultFilename::thermo_data() );
  Antioch::CEAEvaluator<Scalar> thermo_mix( cea_mixture );

// 
  Antioch::PhysicalSet<Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<Scalar> >,
                       Antioch::TransportMixture<Antioch::StatMechThermodynamics<Scalar>,Scalar> >k( tran_mixture );

  Antioch::PhysicalSet<Antioch::BlottnerViscosity<Scalar>,        Antioch::ChemicalMixture<Scalar> > mu( chem_mixture );

  Antioch::PhysicalSet<Antioch::ConstantLewisDiffusivity<Scalar>, Antioch::ChemicalMixture<Scalar> > D( chem_mixture );

// pure species set, all internally set
#ifdef ANTIOCH_HAVE_GSL
  Antioch::PhysicalSet<Antioch::PureSpeciesViscosity<Scalar, Antioch::GSLSpliner>,
                       Antioch::TransportMixture<Antioch::StatMechThermodynamics<Scalar>,Scalar > > ps_mu(tran_mixture);

  Antioch::PhysicalSet<Antioch::MolecularBinaryDiffusion<Scalar, Antioch::GSLSpliner >,
                       Antioch::TransportMixture<Antioch::StatMechThermodynamics<Scalar>,Scalar> > bimol_D( tran_mixture );
#endif

  Antioch::PhysicalSet<Antioch::PureSpeciesThermalConductivity<Antioch::StatMechThermodynamics<Scalar>, Scalar >,
                       Antioch::TransportMixture<Antioch::StatMechThermodynamics<Scalar>,Scalar> > ps_k( tran_mixture );

//Eucken is internally set

  Antioch::read_blottner_data_ascii( mu, Antioch::DefaultFilename::blottner_data() );

  Antioch::build_constant_lewis_diffusivity<Scalar>( D, 1.4);

// non pure species
  Antioch::WilkeMixture< Antioch::PhysicalSet< Antioch::ConstantLewisDiffusivity<Scalar>, Antioch::ChemicalMixture<Scalar> >,
                         Antioch::PhysicalSet< Antioch::BlottnerViscosity<Scalar>,        Antioch::ChemicalMixture<Scalar> >,
                         Antioch::PhysicalSet< Antioch::EuckenThermalConductivity< Antioch::StatMechThermodynamics<Scalar> >,
                                               Antioch::TransportMixture< Antioch::StatMechThermodynamics<Scalar>, Scalar >
                                             >,
                         Antioch::TransportMixture< Antioch::StatMechThermodynamics<Scalar>, Scalar>,
                         Antioch::CEAEvaluator<Scalar>, 
                         Scalar
                       >
        wilke_mixture(D,mu,k, tran_mixture, thermo_mix);

  Antioch::WilkeEvaluator<
                          Antioch::WilkeMixture
                          < 
                           Antioch::PhysicalSet< Antioch::ConstantLewisDiffusivity<Scalar>, Antioch::ChemicalMixture<Scalar> >,
                           Antioch::PhysicalSet< Antioch::BlottnerViscosity<Scalar>,        Antioch::ChemicalMixture<Scalar> >,
                           Antioch::PhysicalSet< Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<Scalar> >, 
                                                 Antioch::TransportMixture<Antioch::StatMechThermodynamics<Scalar>, Scalar >
                                               >,
                           Antioch::TransportMixture< Antioch::StatMechThermodynamics<Scalar>, Scalar>,
                           Antioch::CEAEvaluator<Scalar>, 
                           Scalar
                           >
                         > wilke( wilke_mixture);

// pure species full
#ifdef ANTIOCH_HAVE_GSL
  Antioch::WilkeMixture<
                        Antioch::PhysicalSet<Antioch::MolecularBinaryDiffusion<Scalar, Antioch::GSLSpliner >,
                                             Antioch::TransportMixture<Antioch::StatMechThermodynamics<Scalar>,Scalar> >,
                        Antioch::PhysicalSet<Antioch::PureSpeciesViscosity<Scalar, Antioch::GSLSpliner>,
                                             Antioch::TransportMixture<Antioch::StatMechThermodynamics<Scalar>,Scalar > >,
                        Antioch::PhysicalSet<Antioch::PureSpeciesThermalConductivity<Antioch::StatMechThermodynamics<Scalar>, Scalar >,
                                             Antioch::TransportMixture<Antioch::StatMechThermodynamics<Scalar>,Scalar> >,
                         Antioch::TransportMixture< Antioch::StatMechThermodynamics<Scalar>, Scalar>,
                         Antioch::CEAEvaluator<Scalar>, 
                         Scalar
                       >
        wilke_ps_mixture(bimol_D,ps_mu,ps_k, tran_mixture, thermo_mix);

   Antioch::WilkeEvaluator<
                Antioch::WilkeMixture<
                        Antioch::PhysicalSet<Antioch::MolecularBinaryDiffusion<Scalar, Antioch::GSLSpliner >,
                                             Antioch::TransportMixture<Antioch::StatMechThermodynamics<Scalar>,Scalar> >,
                        Antioch::PhysicalSet<Antioch::PureSpeciesViscosity<Scalar, Antioch::GSLSpliner>,
                                             Antioch::TransportMixture<Antioch::StatMechThermodynamics<Scalar>,Scalar > >,
                        Antioch::PhysicalSet<Antioch::PureSpeciesThermalConductivity<Antioch::StatMechThermodynamics<Scalar>, Scalar >,
                                             Antioch::TransportMixture<Antioch::StatMechThermodynamics<Scalar>,Scalar> >,
                         Antioch::TransportMixture< Antioch::StatMechThermodynamics<Scalar>, Scalar>,
                         Antioch::CEAEvaluator<Scalar>, 
                         Scalar
                               >
                        > wilke_ps_evaluator( wilke_ps_mixture);
#else //only thermal conduction then

  Antioch::WilkeMixture< Antioch::PhysicalSet< Antioch::ConstantLewisDiffusivity<Scalar>, Antioch::ChemicalMixture<Scalar> >,
                         Antioch::PhysicalSet< Antioch::BlottnerViscosity<Scalar>,        Antioch::ChemicalMixture<Scalar> >,
                         Antioch::PhysicalSet< Antioch::PureSpeciesThermalConductivity<Antioch::StatMechThermodynamics<Scalar>, Scalar >,
                                               Antioch::TransportMixture<Antioch::StatMechThermodynamics<Scalar>,Scalar> >,
                         Antioch::TransportMixture< Antioch::StatMechThermodynamics<Scalar>, Scalar>,
                         Antioch::CEAEvaluator<Scalar>, 
                         Scalar
                       >
        wilke_ps_mixture(D,mu,ps_k, tran_mixture, thermo_mix);

  Antioch::WilkeEvaluator<
                          Antioch::WilkeMixture
                          < 
                           Antioch::PhysicalSet< Antioch::ConstantLewisDiffusivity<Scalar>, Antioch::ChemicalMixture<Scalar> >,
                           Antioch::PhysicalSet< Antioch::BlottnerViscosity<Scalar>,        Antioch::ChemicalMixture<Scalar> >,
                           Antioch::PhysicalSet< Antioch::PureSpeciesThermalConductivity<Antioch::StatMechThermodynamics<Scalar>, Scalar >,
                                                 Antioch::TransportMixture<Antioch::StatMechThermodynamics<Scalar>,Scalar> >,
                           Antioch::TransportMixture< Antioch::StatMechThermodynamics<Scalar>, Scalar>,
                           Antioch::CEAEvaluator<Scalar>, 
                           Scalar
                           >
                         > wilke_ps_evaluator( wilke_ps_mixture);
#endif

  int return_flag = 0;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 2;

  // First check phi calculation
  // phi_s = sum_r (chi_r*(1+sqrt(mu_s/mu_r)*(Mr/Ms)^(1/4))^2)/sqrt(8*(1+Ms/Mr))
  {
    std::vector<Scalar> mu( 5 );
    mu[0] = 0.1L;
    mu[1] = 0.2L;
    mu[2] = 0.3L;
    mu[3] = 0.15L;
    mu[4] = 0.25L;

    std::vector<Scalar> chi( 5 );
    chi[0] = 0.1L;
    chi[1] = 0.2L;
    chi[2] = 0.3L;
    chi[3] = 0.15L;
    chi[4] = 0.25L;

    Scalar phi_N_exact = 0.0L;
    unsigned int N_index = 2;
    Scalar M_N = chem_mixture.M(N_index);

    for( unsigned int r = 0; r < 5; r++ )
      {
        Scalar M_r = chem_mixture.M(r);
        Scalar dummy = 1.0L + std::sqrt(mu[N_index]/mu[r])*pow( M_r/M_N, Scalar(0.25L) );
        phi_N_exact += chi[r]*dummy*dummy/std::sqrt(8.0L*( 1.0L + M_N/M_r ) );
      }

    Scalar phi_N;
    phi_N = wilke.compute_phi( mu, chi, N_index );

    return_flag = test_val( phi_N, phi_N_exact, tol, std::string("phi") );
  }
  

  std::vector<Scalar> mass_fractions( 5, 0.2L); 

  // Currently dummy
  //const Scalar mu_exact = ;

  const Scalar T = 1000.0L;
  const Scalar rho = 3.14L;

  Scalar wilke_mu = wilke.mu(T, mass_fractions );
  Scalar wilke_k = wilke.k(T, mass_fractions, rho );
  std::vector<Scalar> wilke_D = Antioch::zero_clone(mass_fractions);
  wilke.D(T, rho, mass_fractions,wilke_D);

  wilke.mu_and_k_and_D(T,rho,mass_fractions,wilke_mu,wilke_k,wilke_D);

  int return_flag_temp = 0;
  //return_flag_temp = test_mu( wilke.mu(T, mass_fractions ), mu_exact, tol );
  if( return_flag_temp != 0 ) return_flag = 1;

  return return_flag;
}

int main()
{
  return ( tester<double>()  ||
           tester<long double>() ||
           tester<float>()
           );
}
