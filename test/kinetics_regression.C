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
#include <string>
#include <vector>

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/chemical_species.h"
#include "antioch/chemical_mixture.h"
#include "antioch/reaction_set.h"
#include "antioch/read_reaction_set_data_xml.h"
#include "antioch/cea_thermo.h"
#include "antioch/kinetics.h"

int main(int argc, char* argv[])
{
  // Check command line count.
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify reaction set XML input file." << std::endl;
      antioch_error();
    }

  std::vector<std::string> species_str_list;
  const unsigned int n_species = 5;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "O2" );
  species_str_list.push_back( "N" );
  species_str_list.push_back( "O" );
  species_str_list.push_back( "NO" );

  Antioch::ChemicalMixture<double> chem_mixture( species_str_list );
  Antioch::ReactionSet<double> reaction_set( chem_mixture );
  Antioch::CEAThermodynamics<double> thermo( chem_mixture );

  std::string input(argv[1]);
  Antioch::read_reaction_set_data_xml<double>( input, true, reaction_set );

  const double T = 1500.0;
  const double P = 1.0e5;

  // Mass fractions
  std::vector<double> Y(n_species,0.2);

  const double R_mix = chem_mixture.R(Y);

  const double rho = P/(R_mix*T);

  std::vector<double> molar_densities(n_species,0.0);
  chem_mixture.molar_densities(rho,Y,molar_densities);

  std::vector<double> h_RT_minus_s_R(n_species);
  thermo.h_RT_minus_s_R(T,h_RT_minus_s_R);

  Antioch::Kinetics<double> kinetics( reaction_set );
  std::vector<double> omega_dot(n_species);
  
  kinetics.compute_mass_sources( T, rho, R_mix, Y, molar_densities, h_RT_minus_s_R, omega_dot );

  for( unsigned int s = 0; s < n_species; s++)
    {
      std::cout << std::scientific << std::setprecision(16)
		<< "omega_dot(" << chem_mixture.chemical_species()[s]->species() << ") = " << omega_dot[s] << std::endl;
    }

  int return_flag = 0;

  std::vector<double> omega_dot_reg(n_species);
  omega_dot_reg[0] =  7.9004530802650654e+04;
  omega_dot_reg[1] = -3.4113853617637843e+05;
  omega_dot_reg[2] = -1.8898881857838202e+05;
  omega_dot_reg[3] =  2.1551399274321867e+05;
  omega_dot_reg[4] =  2.3560883120889112e+05;

  const double tol = 1.0e-14;
  for( unsigned int s = 0; s < n_species; s++)
    {
      const double rel_error = std::fabs( (omega_dot[s] - omega_dot_reg[s])/omega_dot_reg[s]);
      if( rel_error > tol )
	{
	  return_flag = 1;
	}
    }

  if( return_flag == 1 )
    {
      std::cerr << "Error: Mismatch between compute mass source terms and regression values." << std::endl;
      for( unsigned int s = 0; s < n_species; s++)
	{
	  std::cout << std::scientific << std::setprecision(16)
		    << "omega_dot(" << chem_mixture.chemical_species()[s]->species() << ") = " << omega_dot[s]
		    << ", omega_dot_reg(" << chem_mixture.chemical_species()[s]->species() << ") = " << omega_dot_reg[s]
		    << std::endl;
	}
    }
 
  return return_flag;
}
