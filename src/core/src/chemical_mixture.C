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

// This class
#include "antioch/chemical_mixture.h"

// Antioch
#include "antioch/antioch_numeric_type_instantiate_macro.h"
#include "antioch/ascii_parser.h"
#include "antioch/species_parsing.h"

namespace Antioch
{
  template<typename CoeffType>
  ChemicalMixture<CoeffType>::ChemicalMixture(const std::string & filename, const bool verbose,
                                              const std::string & species_data,
                                              const std::string & vibration_data,
                                              const std::string & electronic_data)
  {
    ASCIIParser<CoeffType> parser(filename,verbose);

    read_chemical_species_composition<CoeffType>(static_cast<ParserBase<CoeffType> *> (&parser), *this);

    parser.change_file(species_data);
    this->read_species_characteristics(&parser,species_data,vibration_data,electronic_data);

    return;
  }

  template<typename CoeffType>
  ChemicalMixture<CoeffType>::ChemicalMixture(ParserBase<CoeffType> * parser,
                                              const std::string & species_data,
                                              const std::string & vibration_data,
                                              const std::string & electronic_data)
  {
    read_chemical_species_composition<CoeffType>(parser, *this);

    this->read_species_characteristics(parser,species_data,vibration_data,electronic_data);

    return;
  }

  template<typename CoeffType>
  ChemicalMixture<CoeffType>::ChemicalMixture(const std::vector<std::string>& species_list,
                                              const bool verbose,
                                              const std::string & species_data,
                                              const std::string & vibration_data,
                                              const std::string & electronic_data)
  {
    this->initialize_species(species_list);

    Antioch::ASCIIParser<CoeffType> parser(species_data,verbose);

    this->read_species_characteristics(static_cast<ParserBase<CoeffType> *>(&parser),species_data,vibration_data,electronic_data);

    return;
  }

  template<typename CoeffType>
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
  void ChemicalMixture<CoeffType>::read_species_characteristics(ParserBase<CoeffType> * parser,
                                                                const std::string & /*species_data*/,
                                                                const std::string & vibration_data,
                                                                const std::string & electronic_data)
  {
    // species file is already in parser object
    this->read_species_mandatory_characteristics(parser);

    //... and any vibrational data
    parser->change_file(vibration_data);
    this->read_species_vibrational_characteristics(parser);

    //... and any electronic data
    parser->change_file(electronic_data);
    this->read_species_electronic_characteristics(parser);
  }

  template<typename CoeffType>
  void ChemicalMixture<CoeffType>::initialize_species( const std::vector<std::string>& species_list )
  {
    _chemical_species.resize( species_list.size(), NULL );
    // Build up name map for all possible species
    this->init_species_name_map(species_list);

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
      }
  }

  template<typename CoeffType>
  void ChemicalMixture<CoeffType>::init_species_name_map(const std::vector<std::string> & species_list)
  {
    _species_name_map.clear();
    for(unsigned int s = 0; s < species_list.size(); s++)
      {
        _species_name_map[species_list[s]] = s;
      }
  }

  template<typename CoeffType>
  void ChemicalMixture<CoeffType>::build_inverse_name_map()
  {
    for( std::map<std::string,Species>::const_iterator it = _species_name_map.begin();
	 it != _species_name_map.end(); ++it )
      {
	_species_inv_name_map.insert( std::make_pair( it->second, it->first ) );
      }

    return;
  }

  template <typename CoeffType>
  void ChemicalMixture<CoeffType>::read_species_mandatory_characteristics(ParserBase<CoeffType> * parser)
  {
    read_species_data<CoeffType>(parser, *this);
  }

  template <typename CoeffType>
  void ChemicalMixture<CoeffType>::read_species_vibrational_characteristics(ParserBase<CoeffType> * parser)
  {
    read_species_vibrational_data<CoeffType>(parser, *this);
  }

  template <typename CoeffType>
  void ChemicalMixture<CoeffType>::read_species_electronic_characteristics(ParserBase<CoeffType> * parser)
  {
    read_species_electronic_data<CoeffType>(parser, *this);
  }

  template<typename CoeffType>
  void ChemicalMixture<CoeffType>::add_species( const unsigned int index,
					        const std::string& name,
					        CoeffType mol_wght,
					        CoeffType h_form,
					        CoeffType n_tr_dofs, int charge)
  {
    _chemical_species[index] =
      new ChemicalSpecies<CoeffType>(name, mol_wght, h_form, n_tr_dofs, charge);
  }

  template<typename CoeffType>
  void ChemicalMixture<CoeffType>::add_species_vibrational_data( const unsigned int index,
                                                                 const CoeffType theta_v,
                                                                 const unsigned int ndg_v )
  {
    (_chemical_species[index])->add_vibrational_data(theta_v, ndg_v);
  }

  template<typename CoeffType>
  void ChemicalMixture<CoeffType>::add_species_electronic_data( const unsigned int index,
                                                                const CoeffType theta_e,
                                                                const unsigned int ndg_e )
  {
    (_chemical_species[index])->add_electronic_data(theta_e, ndg_e);
  }

} // end namespace Antioch

// Instantiate
ANTIOCH_NUMERIC_TYPE_CLASS_INSTANTIATE(Antioch::ChemicalMixture);
