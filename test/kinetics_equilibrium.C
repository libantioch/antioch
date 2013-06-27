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

// C++
#include <limits>
#include <string>
#include <vector>

// Antioch
#include "antioch/eigen_utils_decl.h"
#include "antioch/vector_utils.h"

#include "antioch/antioch_asserts.h"
#include "antioch/equilibrium_evaluator.h"
#include "antioch/data_equilibrium.h"
#include "antioch/chemical_species.h"
#include "antioch/chemical_mixture.h"
#include "antioch/reaction_set.h"
#include "antioch/read_reaction_set_data_xml.h"
#include "antioch/cea_thermo.h"
#include "antioch/kinetics_evaluator.h"

template <typename Scalar>
int tester(const std::string& input_name)
{
  using std::abs;

  std::vector<std::string> species_str_list;
  const unsigned int n_species = 5;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "O2" );
  species_str_list.push_back( "NO" );
  species_str_list.push_back( "N" );
  species_str_list.push_back( "O" );

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );
  Antioch::ReactionSet<Scalar> reaction_set( chem_mixture );

  Antioch::read_reaction_set_data_xml<Scalar>( input_name, true, reaction_set );

  const Scalar T = 1500.0;
  const Scalar P = 1.0e5;

  Antioch::KineticsEvaluator<Scalar> kinetics( reaction_set, 0 );


  Antioch::DataEquilibrium<Scalar> equil(T, P,reaction_set);

  std::vector<Scalar> first;
  first.push_back(0.2);
  first.push_back(0.2);
  first.push_back(0.2);
  first.push_back(0.2);
  first.push_back(0.2);

  Antioch::EquilibriumEvaluator<Scalar> eq_solver(equil,kinetics);

  eq_solver.first_guess_molar_fraction(first);
  eq_solver.equilibrium();

  // solution testing
  std::vector<Scalar> Y_eq = eq_solver.mass_fraction_equilibrium();

  const Scalar R_mix = chem_mixture.R(Y_eq);

  std::vector<Scalar> molar_densities = eq_solver.molar_densities_equilibrium();

  Scalar TotDens(0.);
  for(unsigned int nsp = 0; nsp < n_species; nsp++)
  {
     TotDens += molar_densities[nsp];
  }

  const Scalar Peq = eq_solver.Peq();

  const Scalar rho = Peq/(R_mix*T);

  std::vector<Scalar> h_RT_minus_s_R(n_species);
  std::vector<Scalar> dh_RT_minus_s_R_dT(n_species);

  Antioch::CEAThermodynamics<Scalar> thermo( chem_mixture );
  typedef typename Antioch::CEAThermodynamics<Scalar>::template Cache<Scalar> Cache;
  thermo.h_RT_minus_s_R(Cache(T),h_RT_minus_s_R);
  thermo.dh_RT_minus_s_R_dT(Cache(T),dh_RT_minus_s_R_dT);

  std::vector<Scalar> omega_dot(n_species);
  
  kinetics.compute_mass_sources( T, rho, R_mix, Y_eq, molar_densities, h_RT_minus_s_R, omega_dot );

  int return_flag = 0;

  if(eq_solver.success())
  {

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;

  Scalar sum_dot(0.L);
  for( unsigned int s = 0; s < n_species; s++)
    {
      sum_dot += (omega_dot[s] > 0.)?omega_dot[s]:-omega_dot[s];
    }

   if(sum_dot > tol)
   {
      return_flag = 1;
       std::cout << "tolerance is " << tol << " and it is " << sum_dot << std::endl;
   }

   }else
   {
     return_flag = 1;
     std::cout << "equilibrium failed:\n\tsum " << eq_solver.residual()<< ", threshold: " << eq_solver.conv_threshold() << std::endl;
   }

  if(return_flag == 1)
  {
    for( unsigned int s = 0; s < n_species; s++)
      {
        std::cout << std::scientific << std::setprecision(16)
  		  << "omega_dot(" << chem_mixture.chemical_species()[s]->species() << ") = "
		  << omega_dot[s] << "\t"
  		  << "mass_fraction(" << chem_mixture.chemical_species()[s]->species() << ") = "
		  << Y_eq[s] << std::endl;
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

  return (tester<float>(std::string(argv[1])) ||
          tester<double>(std::string(argv[1])) /*||
          tester<long double>(std::string(argv[1])) || */
          );
}
