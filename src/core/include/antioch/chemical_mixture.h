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

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/chemical_species.h"
#include "antioch/metaprogramming.h"
#include "antioch/species_ascii_parsing.h"
#include "antioch/species_enum.h"

// C++
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <iomanip>

namespace Antioch
{
  //! Class storing chemical mixture properties
  /*!
    This class manages the list of ChemicalSpecies for a requested set
    of species from input.
    \todo This should probably be a singleton class, but being lazy for now.
  */
  template<typename CoeffType=double>
  class ChemicalMixture
  {
  public:
    
    ChemicalMixture( const std::vector<std::string>& species_list );
    ~ChemicalMixture();

    //! Returns the number of species in this mixture.
    unsigned int n_species() const;

    void add_species( const unsigned int index, const std::string& name,
		      CoeffType mol_wght, CoeffType h_form,
		      CoeffType n_tr_dofs, int charge );

    void add_species_vibrational_data( const unsigned int index,
                                       const CoeffType theta_v,
                                       const unsigned int ndg_v );

    void add_species_electronic_data( const unsigned int index,
                                      const CoeffType theta_e,
                                      const unsigned int ndg_e );

    const std::vector<ChemicalSpecies<CoeffType>*>& chemical_species() const;

    const std::vector<Species>& species_list() const;

    const std::map<Species,unsigned int>& species_list_map() const;

    const std::map<std::string,unsigned int>& active_species_name_map() const;

    const std::map<std::string,Species>& species_name_map() const;

    const std::map<Species,std::string>& species_inverse_name_map() const;

    //! Gas constant for species s in [J/kg-K]
    CoeffType R( const unsigned int s ) const;

    //! Gas constant for mixture in [J/kg-K]
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type 
    >::type
    R( const VectorStateType& mass_fractions ) const;
    
    //! Molecular weight (molar mass) for species s in [g/mol] or [kg/kmol]
    CoeffType M( const unsigned int s ) const;

    //! Molecular weight (molar mass) for mixture in [g/mol] or [kg/kmol]
    /*!
      \f$ \frac{1}{M} = \sum_s \frac{w_s}{M_s}\f$ where
      \f$ w_s \f$ is the mass fraction of species \f$ s \f$ and
      \f$ M_s \f$ is the molecular weight (molar mass) of species \f$ s \f$
    */
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type 
    >::type
    M( const VectorStateType& mass_fractions ) const;

    //! Species mole fraction
    /*! 
      Given mixture molar mass M and mass fraction for species,
      compute species mole fraction using the relationship
      \f$ w_i = x_i \frac{M_i}{M} \f$ 
    */
    template<typename StateType>
    ANTIOCH_AUTO(StateType) 
    X( const unsigned int species, const StateType& M,
       const StateType& mass_fraction ) const
    ANTIOCH_AUTOFUNC(StateType, mass_fraction*M/this->M(species))

    //! All species mole fractions
    /*!
      The output argument mole_fractions should already be properly
      sized to hold the output.
    */
    template<typename VectorStateType>
    void X( typename Antioch::value_type<VectorStateType>::type M,
	    const VectorStateType& mass_fractions, 
	    VectorStateType& mole_fractions ) const;

    //! Species molar density
    /*! 
      Given total density rho and mass fraction for species,
      compute species moles per unit volume
    */
    template<typename StateType>
    ANTIOCH_AUTO(StateType)
    molar_density( const unsigned int species,
		   const StateType& rho,
		   const StateType& mass_fraction ) const
    ANTIOCH_AUTOFUNC(StateType, rho*mass_fraction/this->M(species))

    //! Species molar densities
    /*! 
      Given total density rho and mass fractions for all species,
      compute moles per unit volume for all species
    */
    template<typename StateType, typename VectorStateType>
    void molar_densities( const StateType& rho,
			  const VectorStateType& mass_fractions,
			  VectorStateType& molar_densities ) const;

  protected:

    void init_species_name_map();
    void build_inverse_name_map();
    void read_species_data();
    void read_species_data( std::istream& in );

    std::vector<Species> _species_list;
    std::map<Species,unsigned int> _species_list_map;
    std::map<std::string,unsigned int> _active_species_name_map;
    std::vector<ChemicalSpecies<CoeffType>*> _chemical_species;
    std::map<std::string,Species> _species_name_map;
    std::map<Species,std::string> _species_inv_name_map;
    
  private:
    ChemicalMixture();

  };


  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType>
  inline
  unsigned int ChemicalMixture<CoeffType>::n_species() const
  {
    return _species_list.size();
  }

  template<typename CoeffType>
  inline
  const std::vector<Species>& ChemicalMixture<CoeffType>::species_list() const
  { 
    return _species_list;
  }

  template<typename CoeffType>
  inline
  const std::map<Species,unsigned int>& ChemicalMixture<CoeffType>::species_list_map() const
  {
    return _species_list_map;
  }

  template<typename CoeffType>
  inline
  const std::map<std::string,unsigned int>& ChemicalMixture<CoeffType>::active_species_name_map() const
  {
    return _active_species_name_map;
  }

  template<typename CoeffType>
  inline
  const std::vector<ChemicalSpecies<CoeffType>*>& ChemicalMixture<CoeffType>::chemical_species() const
  {
    return _chemical_species;
  }

  template<typename CoeffType>
  inline
  const std::map<std::string,Species>& ChemicalMixture<CoeffType>::species_name_map() const
  {
    return _species_name_map;
  }

  template<typename CoeffType>
  inline
  const std::map<Species,std::string>& ChemicalMixture<CoeffType>::species_inverse_name_map() const
  {
    return _species_inv_name_map;
  }

  template<typename CoeffType>
  inline
  CoeffType ChemicalMixture<CoeffType>::R( const unsigned int s ) const
  {
    return (_chemical_species[s])->gas_constant();
  }

  template<typename CoeffType>
  inline
  CoeffType ChemicalMixture<CoeffType>::M( const unsigned int s ) const
  {
    return (_chemical_species[s])->molar_mass();
  }

  template<typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void ChemicalMixture<CoeffType>::molar_densities( const StateType& rho,
						    const VectorStateType& mass_fractions,
						    VectorStateType& molar_densities ) const
  {
    antioch_assert_equal_to( mass_fractions.size(), this->n_species() );
    antioch_assert_equal_to( molar_densities.size(), this->n_species() );

    // !\todo Figure out how to make this assert vector-compatible
    // antioch_assert_greater( rho, 0.0 );

    for( unsigned int s = 0; s < this->n_species(); s++ )
      {
	molar_densities[s] = rho*mass_fractions[s]/this->M(s);
      }
    return;
  }


  template<typename CoeffType>
  inline
  ChemicalMixture<CoeffType>::ChemicalMixture( const std::vector<std::string>& species_list )
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

    //... and any vibrational data
    read_species_vibrational_data_ascii_default(*this);

    //... and any electronic data
    read_species_electronic_data_ascii_default(*this);
  }


  template<typename CoeffType>
  inline
  ChemicalMixture<CoeffType>::~ChemicalMixture()
  {
    // Clean up all the ChemicalSpecies we stored
    for( typename std::vector<ChemicalSpecies<CoeffType>* >::iterator it = _chemical_species.begin();
	 it < _chemical_species.end(); ++it )
      {
	delete (*it);
      }

    return;
  }


  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type 
  >::type
  ChemicalMixture<CoeffType>::R( const VectorStateType& mass_fractions ) const
  {
    antioch_assert_equal_to( mass_fractions.size(), _chemical_species.size() );
    antioch_assert_greater( mass_fractions.size(), 0);
    
    typename Antioch::value_type<VectorStateType>::type 
      R = mass_fractions[0]*this->R(0);
    for( unsigned int s = 1; s < mass_fractions.size(); s++ )
      {
	R += mass_fractions[s]*this->R(s);
      }

    return R;
  }


  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type 
  >::type
  ChemicalMixture<CoeffType>::M( const VectorStateType& mass_fractions ) const
  {
    antioch_assert_equal_to( mass_fractions.size(), _chemical_species.size() );

    typename Antioch::value_type<VectorStateType>::type 
      M = mass_fractions[0]/this->M(0);
    for( unsigned int s = 1; s < mass_fractions.size(); s++ )
      {
	M += mass_fractions[s]/(this->M(s));
      }

    return CoeffType(1)/M;
  }


  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  void ChemicalMixture<CoeffType>::X( typename Antioch::value_type<VectorStateType>::type M,
				      const VectorStateType& mass_fractions, 
				      VectorStateType& mole_fractions ) const
  {
    antioch_assert_equal_to( mass_fractions.size(), _chemical_species.size() );

    // Require output size to be correct.  This is more efficient and
    // works around the difficulty of resizing Eigen containers of
    // types which may lack default constructors.
    antioch_assert_equal_to( mass_fractions.size(), mole_fractions.size() );

    for( unsigned int s = 0; s < mass_fractions.size(); s++ )
      {
	mole_fractions[s] = this->X(s, M, mass_fractions[s]);
      }

    return;
  }


  template<typename CoeffType>
  inline
  void ChemicalMixture<CoeffType>::add_species( const unsigned int index,
					        const std::string& name,
					        CoeffType mol_wght,
					        CoeffType h_form,
					        CoeffType n_tr_dofs, int charge)
  {
    _chemical_species[index] =
      new ChemicalSpecies<CoeffType>(name, mol_wght, h_form, n_tr_dofs, charge);

    return;
  }

  template<typename CoeffType>
  inline
  void ChemicalMixture<CoeffType>::add_species_vibrational_data( const unsigned int index,
                                                                 const CoeffType theta_v,
                                                                 const unsigned int ndg_v )
  {
    (_chemical_species[index])->add_vibrational_data(theta_v, ndg_v);
  }

  template<typename CoeffType>
  inline
  void ChemicalMixture<CoeffType>::add_species_electronic_data( const unsigned int index,
                                                                const CoeffType theta_e,
                                                                const unsigned int ndg_e )
  {
    (_chemical_species[index])->add_electronic_data(theta_e, ndg_e);
  }
  
  template<typename CoeffType>
  inline
  void ChemicalMixture<CoeffType>::init_species_name_map()
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


  template<typename CoeffType>
  inline
  void ChemicalMixture<CoeffType>::build_inverse_name_map()
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
