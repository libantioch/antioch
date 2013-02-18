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
#include "antioch/antioch_asserts.h"
#include "antioch/chemical_species.h"
#include "antioch/species_ascii_parsing.h"
#include "antioch/species_enum.h"

namespace Antioch
{
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

    void add_species( const unsigned int index, const std::string& name,
		      NumericType mol_wght, NumericType h_form,
		      NumericType n_tr_dofs, int charge );

    const std::vector<ChemicalSpecies<NumericType>*>& chemical_species() const;

    const std::vector<Species>& species_list() const;

    const std::map<Species,unsigned int>& species_list_map() const;

    const std::map<std::string,unsigned int>& active_species_name_map() const;

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


  template<class NumericType>
  inline
  ChemicalMixture<NumericType>::ChemicalMixture( const std::vector<std::string>& species_list )
    : _chemical_species( species_list.size(), NULL )
  {
    // Build up name map for all possible species
    this->init_species_name_map();
    
    // Build up inverse name map
    this->build_inverse_name_map();

    // Populate species list for requested species
    _species_list.reserve( species_list.size() );
    for( unsigned int s = 0; s < species_list.size(); s++ )
      {
	if( _species_name_map.find( species_list[s] ) == _species_name_map.end() )
	  {
	    std::cerr << "Error in ChemicalMixture: Unknown species " << species_list[s] << std::endl;
	    antioch_error();
	  }

	_species_list.push_back( _species_name_map.find( species_list[s] )->second );
	_species_list_map.insert( std::make_pair( _species_name_map.find( species_list[s] )->second, s ) );
	_active_species_name_map.insert( std::make_pair( species_list[s], s ) );
      }
	 
    // Now read in chemical properties for the requested species and stash
    read_species_data_ascii_default(*this);
  }


  template<class NumericType>
  inline
  ChemicalMixture<NumericType>::~ChemicalMixture()
  {
    // Clean up all the ChemicalSpecies we stored
    for( typename std::vector<ChemicalSpecies<NumericType>* >::iterator it = _chemical_species.begin();
	 it < _chemical_species.end(); ++it )
      {
	delete (*it);
      }

    return;
  }


  template<class NumericType>
  inline
  NumericType ChemicalMixture<NumericType>::R( const std::vector<NumericType>& mass_fractions ) const
  {
    antioch_assert_equal_to( mass_fractions.size(), _chemical_species.size() );
    
    NumericType R = 0.0;
    for( unsigned int s = 0; s < mass_fractions.size(); s++ )
      {
	R += mass_fractions[s]*this->R(s);
      }

    return R;
  }


  template<class NumericType>
  inline
  NumericType ChemicalMixture<NumericType>::M( const std::vector<NumericType>& mass_fractions ) const
  {
    antioch_assert_equal_to( mass_fractions.size(), _chemical_species.size() );

    NumericType M = 0.0;
    for( unsigned int s = 0; s < mass_fractions.size(); s++ )
      {
	M += mass_fractions[s]/(this->M(s));
      }

    return 1.0/M;
  }


  template<class NumericType>
  inline
  void ChemicalMixture<NumericType>::X( NumericType M,
					const std::vector<NumericType>& mass_fractions, 
					std::vector<NumericType>& mole_fractions ) const
  {
    antioch_assert_equal_to( mass_fractions.size(), _chemical_species.size() );

    mole_fractions.resize( mass_fractions.size() );

    for( unsigned int s = 0; s < mass_fractions.size(); s++ )
      {
	mole_fractions[s] = this->X(s, M, mass_fractions[s]);
      }

    return;
  }


  template<class NumericType>
  inline
  void ChemicalMixture<NumericType>::add_species( const unsigned int index,
						  const std::string& name,
						  NumericType mol_wght,
						  NumericType h_form,
						  NumericType n_tr_dofs, int charge)
  {
    _chemical_species[index] =
      new ChemicalSpecies<NumericType>(name, mol_wght, h_form, n_tr_dofs, charge);

    return;
  }


  template<class NumericType>
  inline
  void ChemicalMixture<NumericType>::init_species_name_map()
  {
    _species_name_map["Air"  ] = Air; 
    _species_name_map["CPAir"] = CPAir; 
    _species_name_map["Ar"   ] = Ar;   
    _species_name_map["Ar+"  ] = Arp;  
    _species_name_map["C"    ] = C;    
    _species_name_map["C+"   ] = Cp;   
    _species_name_map["C2"   ] = C2;   
    _species_name_map["C2H"  ] = C2H;  
    _species_name_map["C2H2" ] = C2H2; 
    _species_name_map["C3"   ] = C3;   
    _species_name_map["CF"   ] = CF;   
    _species_name_map["CF2"  ] = CF2;  
    _species_name_map["CF3"  ] = CF3;  
    _species_name_map["CF4"  ] = CF4;  
    _species_name_map["CH"   ] = CH;   
    _species_name_map["CH2"  ] = CH2;  
    _species_name_map["CH3"  ] = CH3;  
    _species_name_map["CH4"  ] = CH4;  
    _species_name_map["Cl"   ] = Cl;   
    _species_name_map["Cl2"  ] = Cl2;  
    _species_name_map["CN"   ] = CN;   
    _species_name_map["CN+"  ] = CNp;  
    _species_name_map["CO"   ] = CO;   
    _species_name_map["CO+"  ] = COp;  
    _species_name_map["CO2"  ] = CO2;  
    _species_name_map["F"    ] = F;    
    _species_name_map["F2"   ] = F2;   
    _species_name_map["H"    ] = H;    
    _species_name_map["H+"   ] = Hp;   
    _species_name_map["H2"   ] = H2;   
    _species_name_map["H2+"  ] = H2p;  
    _species_name_map["H2O"  ] = H2O;
    _species_name_map["H2O2" ] = H2O2;
    _species_name_map["HCl"  ] = HCl;  
    _species_name_map["HCN"  ] = HCN;  
    _species_name_map["He"   ] = He;   
    _species_name_map["He+"  ] = Hep;
    _species_name_map["HO2"  ] = HO2;
    _species_name_map["N"    ] = N;    
    _species_name_map["N+"   ] = Np;   
    _species_name_map["N2"   ] = N2;   
    _species_name_map["CPN2" ] = CPN2;   
    _species_name_map["N2+"  ] = N2p;  
    _species_name_map["Ne"   ] = Ne;   
    _species_name_map["NCO"  ] = NCO;  
    _species_name_map["NH"   ] = NH;   
    _species_name_map["NH+"  ] = NHp;  
    _species_name_map["NH2"  ] = NH2;  
    _species_name_map["NH3"  ] = NH3;  
    _species_name_map["NO"   ] = NO;   
    _species_name_map["NO+"  ] = NOp;  
    _species_name_map["NO2"  ] = NO2;  
    _species_name_map["O"    ] = O;    
    _species_name_map["O+"   ] = Op;   
    _species_name_map["O2"   ] = O2;   
    _species_name_map["O2+"  ] = O2p;  
    _species_name_map["OH"   ] = OH;   
    _species_name_map["Si"   ] = Si;   
    _species_name_map["SiO"  ] = SiO;  
    _species_name_map["e"    ] = e;

    return;
  }


  template<class NumericType>
  inline
  void ChemicalMixture<NumericType>::build_inverse_name_map()
  {
    for( std::map<std::string,Species>::const_iterator it = _species_name_map.begin();
	 it != _species_name_map.end(); ++it )
      {
	_species_inv_name_map.insert( std::make_pair( it->second, it->first ) );
      }

    return;
  }



} // end namespace Antioch

#endif //ANTIOCH_CHEMICAL_MIXTURE_H
