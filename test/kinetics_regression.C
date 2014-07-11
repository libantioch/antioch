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
#include "antioch/read_reaction_set_data.h"
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
  species_str_list.push_back( "N" );
  species_str_list.push_back( "O" );
  species_str_list.push_back( "NO" );

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );
  Antioch::ReactionSet<Scalar> reaction_set( chem_mixture );
  Antioch::CEAThermodynamics<Scalar> thermo( chem_mixture );

  Antioch::read_reaction_set_data_xml<Scalar>( input_name, true, reaction_set );

  const Scalar T = 1500.0; // K
  const Scalar P = 1.0e5; // Pa

  // Mass fractions
  std::vector<Scalar> Y(n_species,0.2);

  const Scalar R_mix = chem_mixture.R(Y); // get R_tot in J.kg-1.K-1

  const Scalar rho = P/(R_mix*T); // kg.m-3

  std::vector<Scalar> molar_densities(n_species,0.0);
  chem_mixture.molar_densities(rho,Y,molar_densities);

  std::vector<Scalar> h_RT_minus_s_R(n_species);
  std::vector<Scalar> dh_RT_minus_s_R_dT(n_species);

  typedef typename Antioch::CEAThermodynamics<Scalar>::template Cache<Scalar> Cache;
  thermo.h_RT_minus_s_R(Cache(T),h_RT_minus_s_R);
  thermo.dh_RT_minus_s_R_dT(Cache(T),dh_RT_minus_s_R_dT);

  Antioch::KineticsEvaluator<Scalar> kinetics( reaction_set, 0 );

  std::vector<Scalar> omega_dot(n_species);
  std::vector<Scalar> omega_dot_2(n_species);
  std::vector<Scalar> domega_dot_dT(n_species);

  std::vector<std::vector<Scalar> > domega_dot_drho_s(n_species);
  for( unsigned int s = 0; s < n_species; s++ )
    {
      domega_dot_drho_s[s].resize(n_species);
    }
  
  kinetics.compute_mass_sources( T, molar_densities, h_RT_minus_s_R, omega_dot);

  kinetics.compute_mass_sources_and_derivs( T, molar_densities, h_RT_minus_s_R, dh_RT_minus_s_R_dT,
                                            omega_dot_2, domega_dot_dT, domega_dot_drho_s );

  for( unsigned int s = 0; s < n_species; s++)
    {
      std::cout << std::scientific << std::setprecision(16)
                << "domega_dot_dT(" << chem_mixture.chemical_species()[s]->species() << ") = "
                << domega_dot_dT[s] << std::endl;
    }

  for( unsigned int s = 0; s < n_species; s++)
    {
      for( unsigned int t = 0; t < n_species; t++)
        {
          std::cout << std::scientific << std::setprecision(16)
                    << "domega_dot_drho_s(" << chem_mixture.chemical_species()[s]->species() 
                    << ", " << chem_mixture.chemical_species()[t]->species() << ") = "
                    << domega_dot_drho_s[s][t] << std::endl;
        }
    }

  int return_flag = 0;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 500;
  
  // Regression values for omega_dot
  std::vector<Scalar> omega_dot_reg(n_species);
  omega_dot_reg[0] =  9.1623705357123753e+04;
  omega_dot_reg[1] = -3.3462025680272243e+05;
  omega_dot_reg[2] = -2.1139216712069495e+05;
  omega_dot_reg[3] =  1.9782018625609628e+05;
  omega_dot_reg[4] =  2.5656853231019735e+05;
    
  for( unsigned int s = 0; s < n_species; s++)
    {
      const Scalar rel_error = abs( (omega_dot[s] - omega_dot_reg[s])/omega_dot_reg[s]);
      if( rel_error > tol )
        {
          return_flag = 1;
        }
    }
 
  for( unsigned int s = 0; s < n_species; s++)
    {
      const Scalar rel_error = abs( (omega_dot_2[s] - omega_dot_reg[s])/omega_dot_reg[s]);
      if( rel_error > tol )
        {
          return_flag = 1;
        }
    }

  // Regression values for domega_dot_dT
  std::vector<Scalar> domega_dot_reg_dT(n_species);
  domega_dot_reg_dT[0] =  1.8014990183270937e+02;
  domega_dot_reg_dT[1] = -5.2724437115534380e+02;
  domega_dot_reg_dT[2] = -3.0930094476883017e+02;
  domega_dot_reg_dT[3] =  3.7972747459781005e+02;
  domega_dot_reg_dT[4] =  2.7666793949365456e+02;

  for( unsigned int s = 0; s < n_species; s++)
    {
      const Scalar rel_error = abs( (domega_dot_dT[s] - domega_dot_reg_dT[s])/domega_dot_reg_dT[s]);
      if( rel_error > tol )
        {
          return_flag = 1;
        }
    }

  // Regression values for domega_dot_drho_s
  std::vector<std::vector<Scalar> > domega_dot_reg_drhos(n_species);
  for( unsigned int s = 0; s < n_species; s++)
    {
      domega_dot_reg_drhos[s].resize(n_species);
    }

  domega_dot_reg_drhos[0][0] = 1.9675775188085109e+04;
  domega_dot_reg_drhos[0][1] = 1.7226141262419737e+04;
  domega_dot_reg_drhos[0][2] = 3.2159299284723610e+06;
  domega_dot_reg_drhos[0][3] = 1.4765214711933021e+05;
  domega_dot_reg_drhos[0][4] = 2.3225053279918131e+06;

  domega_dot_reg_drhos[1][0] =  8.8927385505978492e+03;
  domega_dot_reg_drhos[1][1] = -9.9560178070099482e+06;
  domega_dot_reg_drhos[1][2] = -9.8748760140991123e+06;
  domega_dot_reg_drhos[1][3] =  4.6143036700500813e+05;
  domega_dot_reg_drhos[1][4] =  8.3487375168772399e+03;

  domega_dot_reg_drhos[2][0] = -2.2420842426881281e+04;
  domega_dot_reg_drhos[2][1] = -4.3812843857644886e+06;
  domega_dot_reg_drhos[2][2] = -6.8343593463263955e+06;
  domega_dot_reg_drhos[2][3] = -5.4143671040862988e+05;
  domega_dot_reg_drhos[2][4] = -1.2267997668149246e+06;

  domega_dot_reg_drhos[3][0] = -1.2028166578920147e+04;
  domega_dot_reg_drhos[3][1] =  4.9713710400172938e+06;
  domega_dot_reg_drhos[3][2] =  5.7418898143800552e+06;
  domega_dot_reg_drhos[3][3] = -9.1121284934572734e+05;
  domega_dot_reg_drhos[3][4] =  1.2431710353864791e+06;

  domega_dot_reg_drhos[4][0] =  5.8804952671184686e+03;
  domega_dot_reg_drhos[4][1] =  9.3487050114947233e+06;
  domega_dot_reg_drhos[4][2] =  7.7514156175730915e+06;
  domega_dot_reg_drhos[4][3] =  8.4356704563001888e+05;
  domega_dot_reg_drhos[4][4] = -2.3472253340802449e+06;

  for( unsigned int s = 0; s < n_species; s++)
    {
      for( unsigned int t = 0; t < n_species; t++)
        {
          const Scalar rel_error = abs( (domega_dot_drho_s[s][t] - domega_dot_reg_drhos[s][t])/domega_dot_reg_drhos[s][t]);
          if( rel_error > tol )
            {
              return_flag = 1;
            }
        }
    }

  // Print out pretty message if there was a problem.
  if( return_flag == 1 )
    {
      std::cerr << "Error: Mismatch between compute mass source terms and regression values." << std::endl;
      for( unsigned int s = 0; s < n_species; s++)
	{
	  std::cout << std::scientific << std::setprecision(16)
		    << "omega_dot(" << chem_mixture.chemical_species()[s]->species() << ") = " << omega_dot[s]
		    << ", omega_dot_reg(" << chem_mixture.chemical_species()[s]->species() << ") = "
                    << omega_dot_reg[s] << std::endl << std::endl;
	}
      
      for( unsigned int s = 0; s < n_species; s++)
	{
	  std::cout << std::scientific << std::setprecision(16)
		    << "domega_dot_dT(" << chem_mixture.chemical_species()[s]->species() << ") = " << domega_dot_dT[s]
		    << ", domega_dot_reg_dT(" << chem_mixture.chemical_species()[s]->species() << ") = "
                    << domega_dot_reg_dT[s] << std::endl << std::endl;
	}

      for( unsigned int s = 0; s < n_species; s++)
	{
          for( unsigned int t = 0; t < n_species; t++)
            {
              std::cout << std::scientific << std::setprecision(16)
                        << "domega_dot_drho_s(" 
                        << chem_mixture.chemical_species()[s]->species() 
                        << ", " << chem_mixture.chemical_species()[t]->species()
                        << ") = " << domega_dot_drho_s[s][t]
                        << ", domega_dot_reg_drho_s("
                        << chem_mixture.chemical_species()[s]->species()
                        << ", " << chem_mixture.chemical_species()[t]->species()
                        << ") = " << domega_dot_reg_drhos[s][t] << std::endl << std::endl;
            }
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
