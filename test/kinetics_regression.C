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

  const Scalar T = 1500.0;
  const Scalar P = 1.0e5;

  // Mass fractions
  std::vector<Scalar> Y(n_species,0.2);

  const Scalar R_mix = chem_mixture.R(Y);

  const Scalar rho = P/(R_mix*T);

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
  
  kinetics.compute_mass_sources( T, rho, R_mix, Y, molar_densities, h_RT_minus_s_R, omega_dot );

  kinetics.compute_mass_sources_and_derivs( T, rho, R_mix, Y, molar_densities, h_RT_minus_s_R, dh_RT_minus_s_R_dT,
                                            omega_dot_2, domega_dot_dT, domega_dot_drho_s );

  for( unsigned int s = 0; s < n_species; s++)
    {
      std::cout << std::scientific << std::setprecision(16)
		<< "omega_dot(" << chem_mixture.chemical_species()[s]->species() << ") = "
		<< omega_dot[s] << std::endl;
    }

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

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;
  
  // Regression values for omega_dot
  std::vector<Scalar> omega_dot_reg(n_species);
  omega_dot_reg[0] =  9.1627505878744108e+04;
  omega_dot_reg[1] = -3.3462516012726480e+05;
  omega_dot_reg[2] = -2.1139750128059302e+05;
  omega_dot_reg[3] =  1.9782333785216609e+05;
  omega_dot_reg[4] =  2.5657181767694763e+05;
    
  for( unsigned int s = 0; s < n_species; s++)
    {
      const Scalar rel_error = std::fabs( (omega_dot[s] - omega_dot_reg[s])/omega_dot_reg[s]);
      if( rel_error > tol )
        {
          return_flag = 1;
        }
    }
 
  for( unsigned int s = 0; s < n_species; s++)
    {
      const Scalar rel_error = std::fabs( (omega_dot_2[s] - omega_dot_reg[s])/omega_dot_reg[s]);
      if( rel_error > tol )
        {
          return_flag = 1;
        }
    }

  // Regression values for domega_dot_dT
  std::vector<Scalar> domega_dot_reg_dT(n_species);
  domega_dot_reg_dT[0] =  1.8015258435766250e+02;
  domega_dot_reg_dT[1] = -5.2724938565430807e+02;
  domega_dot_reg_dT[2] = -3.0930274177637205e+02;
  domega_dot_reg_dT[3] =  3.7973350053870610e+02;
  domega_dot_reg_dT[4] =  2.7666604253431154e+02;

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

  domega_dot_reg_drhos[0][0] = 1.9677528466713775e+04;
  domega_dot_reg_drhos[0][1] = 1.7227676257862015e+04;
  domega_dot_reg_drhos[0][2] = 3.2160890543181915e+06;
  domega_dot_reg_drhos[0][3] = 1.4766530417926086e+05;
  domega_dot_reg_drhos[0][4] = 2.3225848421202623e+06;

  domega_dot_reg_drhos[1][0] =  8.8931540715994815e+03;
  domega_dot_reg_drhos[1][1] = -9.9561695971970223e+06;
  domega_dot_reg_drhos[1][2] = -9.8750240128648076e+06;
  domega_dot_reg_drhos[1][3] =  4.6145192628627969e+05;
  domega_dot_reg_drhos[1][4] =  8.3491261619680736e+03;

  domega_dot_reg_drhos[2][0] = -2.2422758857735978e+04;
  domega_dot_reg_drhos[2][1] = -4.3813526690025711e+06;
  domega_dot_reg_drhos[2][2] = -6.8345704383742884e+06;
  domega_dot_reg_drhos[2][3] = -5.4147327282847988e+05;
  domega_dot_reg_drhos[2][4] = -1.2268436930952054e+06;

  domega_dot_reg_drhos[3][0] = -1.2028768453121129e+04;
  domega_dot_reg_drhos[3][1] =  4.9714465900642881e+06;
  domega_dot_reg_drhos[3][2] =  5.7419784571182663e+06;
  domega_dot_reg_drhos[3][3] = -9.1126114233336027e+05;
  domega_dot_reg_drhos[3][4] =  1.2432112953400959e+06;

  domega_dot_reg_drhos[4][0] =  5.8808447725438464e+03;
  domega_dot_reg_drhos[4][1] =  9.3488479998774435e+06;
  domega_dot_reg_drhos[4][2] =  7.7515269398026383e+06;
  domega_dot_reg_drhos[4][3] =  8.4361718469629961e+05;
  domega_dot_reg_drhos[4][4] = -2.3473015705271205e+06;

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
