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
#include <cmath>
#include <iomanip>
#include <limits>
#include <vector>

// Antioch
#include "antioch/vector_utils.h"

#include "antioch/physical_constants.h"
#include "antioch/chemical_species.h"
#include "antioch/chemical_mixture.h"

template <typename Scalar>
int test_species( const unsigned int species,
		  const std::vector<Antioch::ChemicalSpecies<Scalar>*>& chemical_species,
		  const std::string& species_name,
		  Scalar molar_mass, Scalar gas_constant, Scalar formation_enthalpy, 
		  Scalar n_tr_dofs, int charge,
                  const Scalar tol )
{

  int return_flag = 0;

  const Antioch::ChemicalSpecies<Scalar>& chem_species = *(chemical_species[species]);

  if( chem_species.species() != species_name )
    {
      std::cerr << "Error: Name mismatch for "<< species_name << std::endl
		<< "name = " << chem_species.species() << std::endl;
      return_flag = 1;
    }

  if( std::abs(chem_species.molar_mass() - molar_mass)/molar_mass > tol )
    {
      std::cerr << "Error: Molar mass mismatch for "<< species_name << std::endl
		<< "molar mass = " << chem_species.molar_mass() << std::endl;
      return_flag = 1;
    }

  if( std::abs(chem_species.gas_constant() - gas_constant)/gas_constant > tol )
    {
      std::cerr << "Error: Gas constant mismatch for "<< species_name << std::endl
		<< "gas constant = " << chem_species.gas_constant() << std::endl;
      return_flag = 1;
    }

  if( std::abs(chem_species.formation_enthalpy() - formation_enthalpy)/formation_enthalpy > tol )
    {
      std::cerr << "Error: Formation enthalpy mismatch for "<< species_name << std::endl
		<< "formation enthalpy = " << chem_species.formation_enthalpy() << std::endl;
      return_flag = 1;
    }

  if( chem_species.n_tr_dofs() != n_tr_dofs )
    {
      std::cerr << "Error: Number translational DoFs mismatch for "<< species_name << std::endl
		<< "n_tr_dofs = " << chem_species.n_tr_dofs() << std::endl;
      return_flag = 1;
    }

  if( chem_species.charge() != charge )
    {
      std::cerr << "Error: Charge mismatch for "<< species_name << std::endl
		<< "charge = " << chem_species.charge() << std::endl;
      return_flag = 1;
    }

  return return_flag;
}


template <typename Scalar>
int tester()
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

  const Scalar Mm_N = 14.008e-3L;
  const Scalar Mm_O = 16.000e-3L;
  const Scalar Mm_N2 = 2.L * Mm_N;
  const Scalar Mm_O2 = 2.L * Mm_O;
  const Scalar Mm_NO = Mm_N + Mm_O;

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );
  // testing default
  Antioch::ChemicalMixture<Scalar> default_mixture;
  const unsigned int iN2 = default_mixture.species_name_map().at("N2");
  const unsigned int iO2 = default_mixture.species_name_map().at("O2");
  const unsigned int iNO = default_mixture.species_name_map().at("NO");
  const unsigned int iO  = default_mixture.species_name_map().at("O");
  const unsigned int iN  = default_mixture.species_name_map().at("N");

  const std::map<std::string,Antioch::Species>& species_name_map = chem_mixture.species_name_map();
  const std::map<Antioch::Species,std::string>& species_inverse_name_map = chem_mixture.species_inverse_name_map();
  const std::vector<Antioch::ChemicalSpecies<Scalar>*>& chemical_species = chem_mixture.chemical_species();
  const std::vector<Antioch::Species> species_list = chem_mixture.species_list();
  const std::vector<Antioch::ChemicalSpecies<Scalar>*>& default_species = default_mixture.chemical_species();

  int return_flag = 0;
  const Scalar tol = (std::numeric_limits<Scalar>::epsilon() * 10 < 5e-17)?5e-17:
                                                                           std::numeric_limits<Scalar>::epsilon() * 10;

  // Check name map consistency
  for( unsigned int i = 0; i < n_species; i++ )
    {
      if( species_name_map.find( species_str_list[i] )->second != species_list[i] )
	{
	  std::cerr << "Error: species name map and species list ordering mismatch" << std::endl
		    << "species_name_map = " << species_name_map.find( species_str_list[i] )->second
		    << ", species_list = " << species_list[i] << std::endl;
	  return_flag = 1;
	}
    }

  // Check inverse name map consistency
  for( unsigned int i = 0; i < n_species; i++ )
    {
      if( species_inverse_name_map.find( species_list[i] )->second != species_str_list[i] )
	{
	  std::cerr << "Error: species inverse name map and species list ordering mismatch" << std::endl
		    << "species_inverse_name_map = " << species_inverse_name_map.find( species_list[i] )->second
		    << ", species_str_list = " << species_str_list[i] << std::endl;
	  return_flag = 1;
	}
    }

  // Check N2 properties
  {
    unsigned int index = 0;
    Scalar molar_mass = Mm_N2;
    if( std::abs(molar_mass - chem_mixture.M(index))/molar_mass > tol ||
        std::abs(molar_mass - default_mixture.M(iN2))/molar_mass > tol)
      {
	std::cerr << "Error: Molar mass inconsistency in mixture" << std::endl
		  << "molar mass = " << chem_mixture.M(index) << std::endl
		  << "molar mass (default mixture) = " << default_mixture.M(iN2) << std::endl;
	return_flag = 1;
      }
    return_flag = return_flag || 
                   test_species( index, chemical_species, "N2", molar_mass, 
                                Scalar(Antioch::Constants::R_universal<Scalar>()/molar_mass), 
                                Scalar(0.0), Scalar(2.5), Scalar(0), tol) ||
                   test_species( iN2, default_species, "N2", molar_mass, 
                                Scalar(Antioch::Constants::R_universal<Scalar>()/molar_mass), 
                                Scalar(0.0), Scalar(2.5), Scalar(0), tol);
  }
  
  // Check O2 properties
  {
    unsigned int index = 1;
    Scalar molar_mass = Mm_O2;
    if( std::abs(molar_mass - chem_mixture.M(index)) > tol ||
        std::abs(molar_mass - default_mixture.M(iO2))/molar_mass > tol)
      {
	std::cerr << "Error: Molar mass inconsistency in mixture" << std::endl
		  << "molar mass = " << chem_mixture.M(index) << std::endl
		  << "molar mass (default mixture) = " << default_mixture.M(iO2) << std::endl;
	return_flag = 1;
      }
    return_flag = return_flag || 
                        test_species( index, chemical_species, "O2", molar_mass, 
                                Scalar(Antioch::Constants::R_universal<Scalar>()/molar_mass),
                                Scalar(0.0), Scalar(2.5), Scalar(0), tol) ||
                        test_species( iO2, default_species, "O2", molar_mass, 
                                Scalar(Antioch::Constants::R_universal<Scalar>()/molar_mass),
                                Scalar(0.0), Scalar(2.5), Scalar(0), tol);
  }

  // Check N properties
  {
    unsigned int index = 2;
    Scalar molar_mass = Mm_N;
    if( std::abs(molar_mass - chem_mixture.M(index)) > tol ||
        std::abs(molar_mass - default_mixture.M(iN))/molar_mass > tol)
      {
	std::cerr << "Error: Molar mass inconsistency in mixture" << std::endl
		  << "molar mass = " << chem_mixture.M(index) << std::endl
		  << "molar mass (default mixture) = " << default_mixture.M(iN) << std::endl;
	return_flag = 1;
      }
    return_flag = return_flag || 
                        test_species( index, chemical_species, "N", molar_mass,
                                Scalar(Antioch::Constants::R_universal<Scalar>()/molar_mass), 
                                Scalar(3.3621610000e7), Scalar(1.5), Scalar(0), tol) ||
                        test_species( iN, default_species, "N", molar_mass,
                                Scalar(Antioch::Constants::R_universal<Scalar>()/molar_mass), 
                                Scalar(3.3621610000e7), Scalar(1.5), Scalar(0), tol);
  }

  // Check O properties
  {
    unsigned int index = 3;
    Scalar molar_mass = Mm_O;
    if( std::abs(molar_mass - chem_mixture.M(index)) > tol ||
        std::abs(molar_mass - default_mixture.M(iO)) > tol )
      {
	std::cerr << "Error: Molar mass inconsistency in mixture" << std::endl
		  << "molar mass = " << chem_mixture.M(index) << std::endl
		  << "molar mass (default mixture) = " << default_mixture.M(iO) << std::endl;
	return_flag = 1;
      }
    return_flag = return_flag || 
                        test_species( index, chemical_species, "O", molar_mass,
                                Scalar(Antioch::Constants::R_universal<Scalar>()/molar_mass), 
                                Scalar(1.5420000000e7), Scalar(1.5), Scalar(0), tol) ||
                        test_species( iO, default_species, "O", molar_mass,
                                Scalar(Antioch::Constants::R_universal<Scalar>()/molar_mass), 
                                Scalar(1.5420000000e7), Scalar(1.5), Scalar(0), tol);
  }

  // Check NO properties
  {
    unsigned int index = 4;
    Scalar molar_mass = Mm_NO;
    if( std::abs(molar_mass - chem_mixture.M(index)) > tol ||
        std::abs(molar_mass - default_mixture.M(iNO)) > tol )
      {
	std::cerr << "Error: Molar mass inconsistency in mixture" << std::endl
		  << "molar mass = " << chem_mixture.M(index) << std::endl
		  << "molar mass (default mixture) = " << default_mixture.M(iNO) << std::endl;
	return_flag = 1;
      }
    return_flag = return_flag || 
                        test_species( index, chemical_species, "NO", molar_mass,
                                Scalar(Antioch::Constants::R_universal<Scalar>()/molar_mass), 
                                Scalar(2.9961230000e6), Scalar(2.5), Scalar(0), tol) ||
                        test_species( iNO, default_species, "NO", molar_mass,
                                Scalar(Antioch::Constants::R_universal<Scalar>()/molar_mass), 
                                Scalar(2.9961230000e6), Scalar(2.5), Scalar(0), tol);
  }

  std::vector<Scalar> mass_fractions( 5, 0.2L );

  Scalar R_exact = Antioch::Constants::R_universal<Scalar>()*( 0.2L/Mm_N2 + 0.2L/Mm_O2 + 0.2L/Mm_N + 0.2L/Mm_O + 0.2L/Mm_NO );
  Scalar M_exact = 1.0L/( 0.2L*( 1.0L/Mm_N2 + 1.0L/Mm_O2 + 1.0L/Mm_N + 1.0L/Mm_O + 1.0L/Mm_NO) );
  
  std::vector<Scalar> X_exact(5, 0.0L);
  X_exact[0] = 0.2L*M_exact/Mm_N2;
  X_exact[1] = 0.2L*M_exact/Mm_O2;
  X_exact[2] = 0.2L*M_exact/Mm_N;
  X_exact[3] = 0.2L*M_exact/Mm_O;
  X_exact[4] = 0.2L*M_exact/Mm_NO;

  if( abs( (chem_mixture.R(mass_fractions) - R_exact)/R_exact) > tol )
    {
      std::cerr << "Error: Mismatch in mixture gas constant." << std::endl
		<< std::setprecision(16) << std::scientific
		<< "R       = " << chem_mixture.R(mass_fractions) << std::endl
		<< "R_exact = " << R_exact <<  std::endl;
      return_flag = 1;
    }

  if( abs( (chem_mixture.M(mass_fractions) - M_exact)/M_exact ) > tol )
    {
      std::cerr << "Error: Mismatch in mixture molar mass." << std::endl
		<< std::setprecision(16) << std::scientific
		<< "M       = " << chem_mixture.M(mass_fractions) << std::endl
		<< "M_exact = " << M_exact << std::endl;
      return_flag = 1;
    }
  
  std::vector<Scalar> X(5, 0);
  chem_mixture.X( chem_mixture.M(mass_fractions), mass_fractions, X );
  for( unsigned int s = 0; s < 5; s++ )
    {
      if( abs( (X[s] - X_exact[s])/X_exact[s]) > tol )
	{
	  std::cerr << "Error: Mismatch in mole fraction for species " << s << std::endl
		    << std::setprecision(16) << std::scientific
		    << "X       = " << X[s] << std::endl
		    << "X_exact = " << X_exact[s] << std::endl;
	  return_flag = 1;
	}
    }

  return return_flag;
}


int main()
{
  return (tester<double>() ||
          tester<long double>() ||
          tester<float>());
}
