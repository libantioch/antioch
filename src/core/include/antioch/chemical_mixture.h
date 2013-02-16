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

#ifndef ANTIOCH_CHEMICAL_MIXTURE_H
#define ANTIOCH_CHEMICAL_MIXTURE_H

// C++
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <iomanip>

// Antioch
#include "antioch/species_enum.h"
#include "antioch/antioch_asserts.h"

namespace Antioch
{
  // Forward declarations
  template<class NumericType>
  class ChemicalSpecies;

  //! Class storing chemical mixture properties
  /*!
    This class manages the list of ChemicalSpecies for a requested set
    of species from input.
    \todo This should probably be a singleton class, but being lazy for now.
  */
  template<class NumericType>
  class ChemicalMixture
  {
  public:
    
    ChemicalMixture( const std::vector<std::string>& species_list );
    ~ChemicalMixture();

    //! Returns the number of species in this mixture.
    unsigned int n_species() const;

    const std::vector<Species>& species_list() const;

    const std::map<Species,unsigned int>& species_list_map() const;

    const std::map<std::string,unsigned int>& active_species_name_map() const;

    const std::vector<ChemicalSpecies<NumericType>*>& chemical_species() const;

    const std::map<std::string,Species>& species_name_map() const;

    const std::map<Species,std::string>& species_inverse_name_map() const;

    //! Gas constant for species s in [J/kg-K]
    NumericType R( const unsigned int s ) const;

    //! Gas constant for mixture in [J/kg-K]
    NumericType R( const std::vector<NumericType>& mass_fractions ) const;
    
    //! Molecular weight (molar mass) for species s in [g/mol] or [kg/kmol]
    NumericType M( const unsigned int s ) const;

    //! Molecular weight (molar mass) for mixture in [g/mol] or [kg/kmol]
    /*!
      \f$ \frac{1}{M} = \sum_s \frac{w_s}{M_s}\f$ where
      \f$ w_s \f$ is the mass fraction of species \f$ s \f$ and
      \f$ M_s \f$ is the molecular weight (molar mass) of species \f$ s \f$
    */
    NumericType M( const std::vector<NumericType>& mass_fractions ) const;

    //! Species mole fraction
    /*! 
      Given mixture molar mass M and mass fraction for species,
      compute species mole fraction using the relationship
      \f$ w_i = x_i \frac{M_i}{M} \f$ 
    */
    NumericType X( const unsigned int species, const NumericType M,
		   const NumericType mass_fraction ) const;

    //! All species mole fractions
    void X( NumericType M, const std::vector<NumericType>& mass_fractions, 
	    std::vector<NumericType>& mole_fractions ) const;

    NumericType molar_density( const unsigned int species,
			       const NumericType rho,
			       const NumericType mass_fraction ) const;

    void molar_densities( const NumericType rho,
			  const std::vector<NumericType>& mass_fractions,
			  std::vector<NumericType>& molar_densities ) const;

  protected:

    void init_species_name_map();
    void build_inverse_name_map();
    void read_species_data();
    void read_species_data( std::istream& in );

    std::vector<Species> _species_list;
    std::map<Species,unsigned int> _species_list_map;
    std::map<std::string,unsigned int> _active_species_name_map;
    std::vector<ChemicalSpecies<NumericType>*> _chemical_species;
    std::map<std::string,Species> _species_name_map;
    std::map<Species,std::string> _species_inv_name_map;
    
  private:
    ChemicalMixture();

  };

  /* ------------------------- Inline Functions -------------------------*/
  template<class NumericType>
  inline
  unsigned int ChemicalMixture<NumericType>::n_species() const
  {
    return _species_list.size();
  }

  template<class NumericType>
  inline
  const std::vector<Species>& ChemicalMixture<NumericType>::species_list() const
  { 
    return _species_list;
  }

  template<class NumericType>
  inline
  const std::map<Species,unsigned int>& ChemicalMixture<NumericType>::species_list_map() const
  {
    return _species_list_map;
  }

  template<class NumericType>
  inline
  const std::map<std::string,unsigned int>& ChemicalMixture<NumericType>::active_species_name_map() const
  {
    return _active_species_name_map;
  }

  template<class NumericType>
  inline
  const std::vector<ChemicalSpecies<NumericType>*>& ChemicalMixture<NumericType>::chemical_species() const
  {
    return _chemical_species;
  }

  template<class NumericType>
  inline
  const std::map<std::string,Species>& ChemicalMixture<NumericType>::species_name_map() const
  {
    return _species_name_map;
  }

  template<class NumericType>
  inline
  const std::map<Species,std::string>& ChemicalMixture<NumericType>::species_inverse_name_map() const
  {
    return _species_inv_name_map;
  }

  template<class NumericType>
  inline
  NumericType ChemicalMixture<NumericType>::R( const unsigned int s ) const
  {
    return (_chemical_species[s])->gas_constant();
  }

  template<class NumericType>
  inline
  NumericType ChemicalMixture<NumericType>::M( const unsigned int s ) const
  {
    return (_chemical_species[s])->molar_mass();
  }

  template<class NumericType>
  inline
  NumericType ChemicalMixture<NumericType>::X( const unsigned int species,
					       const NumericType M,
					       const NumericType mass_fraction ) const
  {
    return mass_fraction*M/this->M(species);
  }

  template<class NumericType>
  inline
  NumericType ChemicalMixture<NumericType>::molar_density( const unsigned int species,
							   const NumericType rho,
							   const NumericType mass_fraction ) const
  {
    antioch_assert_greater( rho, 0.0 );
    return rho*mass_fraction/this->M(species);
  }

  template<class NumericType>
  inline
  void ChemicalMixture<NumericType>::molar_densities( const NumericType rho,
						      const std::vector<NumericType>& mass_fractions,
						      std::vector<NumericType>& molar_densities ) const
  {
    antioch_assert_equal_to( mass_fractions.size(), this->n_species() );
    antioch_assert_equal_to( molar_densities.size(), this->n_species() );
    antioch_assert_greater( rho, 0.0 );
    for( unsigned int s = 0; s < this->n_species(); s++ )
      {
	molar_densities[s] = rho*mass_fractions[s]/this->M(s);
      }
    return;
  }

} // end namespace Antioch

#endif //ANTIOCH_CHEMICAL_MIXTURE_H
