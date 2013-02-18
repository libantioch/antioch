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
  template<typename CoefType=double>
  class ChemicalMixture
  {
  public:
    
    ChemicalMixture( const std::vector<std::string>& species_list );
    ~ChemicalMixture();

    //! Returns the number of species in this mixture.
    unsigned int n_species() const;

    void add_species( const unsigned int index, const std::string& name,
		      CoefType mol_wght, CoefType h_form,
		      CoefType n_tr_dofs, int charge );

    const std::vector<ChemicalSpecies<CoefType>*>& chemical_species() const;

    const std::vector<Species>& species_list() const;

    const std::map<Species,unsigned int>& species_list_map() const;

    const std::map<std::string,unsigned int>& active_species_name_map() const;

    const std::map<std::string,Species>& species_name_map() const;

    const std::map<Species,std::string>& species_inverse_name_map() const;

    //! Gas constant for species s in [J/kg-K]
    CoefType R( const unsigned int s ) const;

    //! Gas constant for mixture in [J/kg-K]
    template<typename StateType>
    StateType R( const std::vector<StateType>& mass_fractions ) const;
    
    //! Molecular weight (molar mass) for species s in [g/mol] or [kg/kmol]
    CoefType M( const unsigned int s ) const;

    //! Molecular weight (molar mass) for mixture in [g/mol] or [kg/kmol]
    /*!
      \f$ \frac{1}{M} = \sum_s \frac{w_s}{M_s}\f$ where
      \f$ w_s \f$ is the mass fraction of species \f$ s \f$ and
      \f$ M_s \f$ is the molecular weight (molar mass) of species \f$ s \f$
    */
    template<typename StateType>
    StateType M( const std::vector<StateType>& mass_fractions ) const;

    //! Species mole fraction
    /*! 
      Given mixture molar mass M and mass fraction for species,
      compute species mole fraction using the relationship
      \f$ w_i = x_i \frac{M_i}{M} \f$ 
    */
    template<typename StateType>
    StateType X( const unsigned int species, const StateType M,
		 const StateType mass_fraction ) const;

    //! All species mole fractions
    template<typename StateType>
    void X( StateType M, const std::vector<StateType>& mass_fractions, 
	    std::vector<StateType>& mole_fractions ) const;

    //! Species molar density
    /*! 
      Given total density rho and mass fraction for species,
      compute species moles per unit volume
    */
    template<typename StateType>
    StateType molar_density( const unsigned int species,
			     const StateType rho,
			     const StateType mass_fraction ) const;

    //! Species molar densities
    /*! 
      Given total density rho and mass fractions for all species,
      compute moles per unit volume for all species
    */
    template<typename StateType>
    void molar_densities( const StateType rho,
			  const std::vector<StateType>& mass_fractions,
			  std::vector<StateType>& molar_densities ) const;

  protected:

    void init_species_name_map();
    void build_inverse_name_map();
    void read_species_data();
    void read_species_data( std::istream& in );

    std::vector<Species> _species_list;
    std::map<Species,unsigned int> _species_list_map;
    std::map<std::string,unsigned int> _active_species_name_map;
    std::vector<ChemicalSpecies<CoefType>*> _chemical_species;
    std::map<std::string,Species> _species_name_map;
    std::map<Species,std::string> _species_inv_name_map;
    
  private:
    ChemicalMixture();

  };


  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoefType>
  inline
  unsigned int ChemicalMixture<CoefType>::n_species() const
  {
    return _species_list.size();
  }

  template<typename CoefType>
  inline
  const std::vector<Species>& ChemicalMixture<CoefType>::species_list() const
  { 
    return _species_list;
  }

  template<typename CoefType>
  inline
  const std::map<Species,unsigned int>& ChemicalMixture<CoefType>::species_list_map() const
  {
    return _species_list_map;
  }

  template<typename CoefType>
  inline
  const std::map<std::string,unsigned int>& ChemicalMixture<CoefType>::active_species_name_map() const
  {
    return _active_species_name_map;
  }

  template<typename CoefType>
  inline
  const std::vector<ChemicalSpecies<CoefType>*>& ChemicalMixture<CoefType>::chemical_species() const
  {
    return _chemical_species;
  }

  template<typename CoefType>
  inline
  const std::map<std::string,Species>& ChemicalMixture<CoefType>::species_name_map() const
  {
    return _species_name_map;
  }

  template<typename CoefType>
  inline
  const std::map<Species,std::string>& ChemicalMixture<CoefType>::species_inverse_name_map() const
  {
    return _species_inv_name_map;
  }

  template<typename CoefType>
  inline
  CoefType ChemicalMixture<CoefType>::R( const unsigned int s ) const
  {
    return (_chemical_species[s])->gas_constant();
  }

  template<typename CoefType>
  inline
  CoefType ChemicalMixture<CoefType>::M( const unsigned int s ) const
  {
    return (_chemical_species[s])->molar_mass();
  }

  template<typename CoefType>
  template<typename StateType>
  inline
  StateType ChemicalMixture<CoefType>::X( const unsigned int species,
					  const StateType M,
					  const StateType mass_fraction ) const
  {
    return mass_fraction*M/this->M(species);
  }

  template<typename CoefType>
  template<typename StateType>
  inline
  StateType ChemicalMixture<CoefType>::molar_density( const unsigned int species,
						      const StateType rho,
						      const StateType mass_fraction ) const
  {
    antioch_assert_greater( rho, 0.0 );
    return rho*mass_fraction/this->M(species);
  }

  template<typename CoefType>
  template<typename StateType>
  inline
  void ChemicalMixture<CoefType>::molar_densities( const StateType rho,
						   const std::vector<StateType>& mass_fractions,
						   std::vector<StateType>& molar_densities ) const
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


  template<typename CoefType>
  inline
  ChemicalMixture<CoefType>::ChemicalMixture( const std::vector<std::string>& species_list )
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


  template<typename CoefType>
  inline
  ChemicalMixture<CoefType>::~ChemicalMixture()
  {
    // Clean up all the ChemicalSpecies we stored
    for( typename std::vector<ChemicalSpecies<CoefType>* >::iterator it = _chemical_species.begin();
	 it < _chemical_species.end(); ++it )
      {
	delete (*it);
      }

    return;
  }


  template<typename CoefType>
  template<typename StateType>
  inline
  StateType ChemicalMixture<CoefType>::R( const std::vector<StateType>& mass_fractions ) const
  {
    antioch_assert_equal_to( mass_fractions.size(), _chemical_species.size() );
    
    StateType R = 0.0;
    for( unsigned int s = 0; s < mass_fractions.size(); s++ )
      {
	R += mass_fractions[s]*this->R(s);
      }

    return R;
  }


  template<typename CoefType>
  template<typename StateType>
  inline
  StateType ChemicalMixture<CoefType>::M( const std::vector<StateType>& mass_fractions ) const
  {
    antioch_assert_equal_to( mass_fractions.size(), _chemical_species.size() );

    StateType M = 0.0;
    for( unsigned int s = 0; s < mass_fractions.size(); s++ )
      {
	M += mass_fractions[s]/(this->M(s));
      }

    return 1.0/M;
  }


  template<typename CoefType>
  template<typename StateType>
  inline
  void ChemicalMixture<CoefType>::X( StateType M,
				     const std::vector<StateType>& mass_fractions, 
				     std::vector<StateType>& mole_fractions ) const
  {
    antioch_assert_equal_to( mass_fractions.size(), _chemical_species.size() );

    mole_fractions.resize( mass_fractions.size() );

    for( unsigned int s = 0; s < mass_fractions.size(); s++ )
      {
	mole_fractions[s] = this->X(s, M, mass_fractions[s]);
      }

    return;
  }


  template<typename CoefType>
  inline
  void ChemicalMixture<CoefType>::add_species( const unsigned int index,
					       const std::string& name,
					       CoefType mol_wght,
					       CoefType h_form,
					       CoefType n_tr_dofs, int charge)
  {
    _chemical_species[index] =
      new ChemicalSpecies<CoefType>(name, mol_wght, h_form, n_tr_dofs, charge);

    return;
  }


  template<typename CoefType>
  inline
  void ChemicalMixture<CoefType>::init_species_name_map()
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


  template<typename CoefType>
  inline
  void ChemicalMixture<CoefType>::build_inverse_name_map()
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
