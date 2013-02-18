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
#include <sstream>

// This class
#include "antioch/chemical_mixture.h"

// Antioch
#include "antioch/input_utils.h"
#include "antioch/chemical_species.h"
#include "antioch/species_ascii_parsing.h"

namespace Antioch
{
  template<class NumericType>
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
  void ChemicalMixture<NumericType>::build_inverse_name_map()
  {
    for( std::map<std::string,Species>::const_iterator it = _species_name_map.begin();
	 it != _species_name_map.end(); ++it )
      {
	_species_inv_name_map.insert( std::make_pair( it->second, it->first ) );
      }

    return;
  }

  /* ------------------------- Instantiate ------------------------- */
  template class ChemicalMixture<double>;

} // end namespace Antioch
