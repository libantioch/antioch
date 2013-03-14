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
#include <cmath>
#include <iomanip>
#include <limits>

// Antioch
#include "antioch/eigen_utils.h"
#include "antioch/metaphysicl_utils.h"
#include "antioch/valarray_utils.h"
#include "antioch/vector_utils.h"

#include "antioch/chemical_mixture.h"
#include "antioch/chemical_species.h"
#include "antioch/physical_constants.h"


template <typename Scalar>
int test_species( const unsigned int species,
		  const std::vector<Antioch::ChemicalSpecies<Scalar>*>& chemical_species,
		  const std::string& species_name,
		  Scalar molar_mass, Scalar gas_constant, Scalar formation_enthalpy, 
		  Scalar n_tr_dofs, int charge )
{

  int return_flag = 0;

  const Antioch::ChemicalSpecies<Scalar>& chem_species = *(chemical_species[species]);

  if( chem_species.species() != species_name )
    {
      std::cerr << "Error: Name mismatch for "<< species_name << std::endl
		<< "name = " << chem_species.species() << std::endl;
      return_flag = 1;
    }

  if( chem_species.molar_mass() != molar_mass )
    {
      std::cerr << "Error: Molar mass mismatch for "<< species_name << std::endl
		<< "molar mass = " << chem_species.molar_mass() << std::endl;
      return_flag = 1;
    }

  if( chem_species.gas_constant() != gas_constant )
    {
      std::cerr << "Error: Gas constant mismatch for "<< species_name << std::endl
		<< "gas constant = " << chem_species.gas_constant() << std::endl;
      return_flag = 1;
    }

  if( chem_species.formation_enthalpy() != formation_enthalpy )
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


template <typename Scalar, typename PairScalars>
int vectester(const PairScalars& example)
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

  const std::map<std::string,Antioch::Species>& species_name_map = chem_mixture.species_name_map();
  const std::map<Antioch::Species,std::string>& species_inverse_name_map = chem_mixture.species_inverse_name_map();
  const std::vector<Antioch::ChemicalSpecies<Scalar>*>& chemical_species = chem_mixture.chemical_species();
  const std::vector<Antioch::Species> species_list = chem_mixture.species_list();

  int return_flag = 0;

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
    Scalar molar_mass = 28.01600L;
    if( molar_mass != chem_mixture.M(index) )
      {
	std::cerr << "Error: Molar mass inconsistency in mixture" << std::endl
		  << "molar mass = " << chem_mixture.M(index) << std::endl;
	return_flag = 1;
      }
    return_flag = test_species( index, chemical_species, "N2", molar_mass, 
                                Scalar(Antioch::Constants::R_universal/molar_mass), 
                                Scalar(0.0), Scalar(2.5), Scalar(0));
  }
  
  // Check O2 properties
  {
    unsigned int index = 1;
    Scalar molar_mass = 32.00000L;
    if( molar_mass != chem_mixture.M(index) )
      {
	std::cerr << "Error: Molar mass inconsistency in mixture" << std::endl
		  << "molar mass = " << chem_mixture.M(index) << std::endl;
	return_flag = 1;
      }
    return_flag = test_species( index, chemical_species, "O2", molar_mass, 
                                Scalar(Antioch::Constants::R_universal/molar_mass),
                                Scalar(0.0), Scalar(2.5), Scalar(0));
  }

  // Check N properties
  {
    unsigned int index = 2;
    Scalar molar_mass = 14.00800L;
    if( molar_mass != chem_mixture.M(index) )
      {
	std::cerr << "Error: Molar mass inconsistency in mixture" << std::endl
		  << "molar mass = " << chem_mixture.M(index) << std::endl;
	return_flag = 1;
      }
    return_flag = test_species( index, chemical_species, "N", molar_mass,
                                Scalar(Antioch::Constants::R_universal/molar_mass), 
                                Scalar(3.3621610000e7), Scalar(1.5), Scalar(0));
  }

  // Check O properties
  {
    unsigned int index = 3;
    Scalar molar_mass = 16.00000L;
    if( molar_mass != chem_mixture.M(index) )
      {
	std::cerr << "Error: Molar mass inconsistency in mixture" << std::endl
		  << "molar mass = " << chem_mixture.M(index) << std::endl;
	return_flag = 1;
      }
    return_flag = test_species( index, chemical_species, "O", molar_mass,
                                Scalar(Antioch::Constants::R_universal/molar_mass), 
                                Scalar(1.5420000000e7), Scalar(1.5), Scalar(0));
  }

  // Check NO properties
  {
    unsigned int index = 4;
    Scalar molar_mass = 30.00800L;
    if( molar_mass != chem_mixture.M(index) )
      {
	std::cerr << "Error: Molar mass inconsistency in mixture" << std::endl
		  << "molar mass = " << chem_mixture.M(index) << std::endl;
	return_flag = 1;
      }
    return_flag = test_species( index, chemical_species, "NO", molar_mass,
                                Scalar(Antioch::Constants::R_universal/molar_mass), 
                                Scalar(2.9961230000e6), Scalar(2.5), Scalar(0));
  }

  std::vector<PairScalars> mass_fractions( 5, example );
  mass_fractions[0][0] = 0.25;
  mass_fractions[1][0] = 0.25;
  mass_fractions[2][0] = 0.25;
  mass_fractions[3][0] = 0.25;
  mass_fractions[4][0] = 0;
  mass_fractions[0][1] = 0.2;
  mass_fractions[1][1] = 0.2;
  mass_fractions[2][1] = 0.2;
  mass_fractions[3][1] = 0.2;
  mass_fractions[4][1] = 0.2;

  PairScalars R_exact = example;
  PairScalars M_exact = example;
  R_exact[0] = Antioch::Constants::R_universal*( 0.25/28.016 + 0.25/32.0 + 0.25/14.008 + 0.25/16.0);
  R_exact[1] = Antioch::Constants::R_universal*( 0.2/28.016 + 0.2/32.0 + 0.2/14.008 + 0.2/16.0 + 0.2/30.008 );
  M_exact[0] = 1.0/( 0.25*( 1.0/28.016 + 1.0/32.0 + 1.0/14.008 + 1.0/16.0) );
  M_exact[1] = 1.0/( 0.2*( 1.0/28.016 + 1.0/32.0 + 1.0/14.008 + 1.0/16.0 + 1.0/30.008) );
  
  std::vector<PairScalars> X_exact(5, example);
  X_exact[0][0] = 0.25*M_exact[0]/28.016;
  X_exact[1][0] = 0.25*M_exact[0]/32.0;
  X_exact[2][0] = 0.25*M_exact[0]/14.008;
  X_exact[3][0] = 0.25*M_exact[0]/16.0;
  X_exact[4][0] = 0;
  X_exact[0][1] = 0.2*M_exact[1]/28.016;
  X_exact[1][1] = 0.2*M_exact[1]/32.0;
  X_exact[2][1] = 0.2*M_exact[1]/14.008;
  X_exact[3][1] = 0.2*M_exact[1]/16.0;
  X_exact[4][1] = 0.2*M_exact[1]/30.008;

  Scalar tol = std::numeric_limits<Scalar>::epsilon() * 10;
  const PairScalars rel_R_error = 
    std::abs( (chem_mixture.R(mass_fractions) - R_exact)/R_exact);
  if( Antioch::max(rel_R_error) > tol )
    {
      std::cerr << "Error: Mismatch in mixture gas constant." << std::endl
		<< std::setprecision(16) << std::scientific
		<< "R       = " << chem_mixture.R(mass_fractions) << std::endl
		<< "R_exact = " << R_exact <<  std::endl;
      return_flag = 1;
    }

  const PairScalars rel_M_error = 
    std::abs( (chem_mixture.M(mass_fractions) - M_exact)/M_exact);
  if( Antioch::max(rel_M_error) > tol )
    {
      std::cerr << "Error: Mismatch in mixture molar mass." << std::endl
		<< std::setprecision(16) << std::scientific
		<< "M       = " << chem_mixture.M(mass_fractions) << std::endl
		<< "M_exact = " << M_exact << std::endl;
      return_flag = 1;
    }
  
  std::vector<PairScalars> X;
  chem_mixture.X( chem_mixture.M(mass_fractions), mass_fractions, X );
  for( unsigned int s = 0; s < 5; s++ )
    {
      const PairScalars rel_X_error = 
        std::abs( (X[s] - X_exact[s])/X_exact[s]);
      if( Antioch::max(rel_X_error) > tol )
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
  int returnval = 0;

  returnval = returnval ||
    vectester<float, std::valarray<float> >
      (std::valarray<float>(2));
  returnval = returnval ||
    vectester<double, std::valarray<double> >
      (std::valarray<double>(2));
// We're not getting the full long double precision yet?
//  returnval = returnval ||
//    vectester<long double, std::valarray<long double> >
//      (std::valarray<long double>(2));
#ifdef ANTIOCH_HAVE_EIGEN
  returnval = returnval ||
    vectester<float, Eigen::Array2f>
      (Eigen::Array2f());
  returnval = returnval ||
    vectester<double, Eigen::Array2d>
      (Eigen::Array2d());
// We're not getting the full long double precision yet?
//  returnval = returnval ||
//    vectester<long double, Eigen::Array<long double, 2, 1> >
//      (Eigen::Array<long double, 2, 1>());
#endif
#ifdef ANTIOCH_HAVE_METAPHYSICL
  returnval = returnval ||
    vectester<float, MetaPhysicL::NumberArray<2, float> > (0);
  returnval = returnval ||
    vectester<double, MetaPhysicL::NumberArray<2, double> > (0);
//  returnval = returnval ||
//    vectester<long double, MetaPhysicL::NumberArray<2, long double> > (0);
#endif

  return returnval;
}
