//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014-2016 Paul T. Bauman, Benjamin S. Kirk,
//                         Sylvain Plessis, Roy H. Stonger
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

// This class
#include "antioch/transport_mixture.h"

// Antioch
#include "antioch/antioch_numeric_type_instantiate_macro.h"
#include "antioch/ascii_parser.h"
#include "antioch/xml_parser.h"
#include "antioch/chemkin_parser.h"

namespace Antioch
{
  template<typename CoeffType>
  TransportMixture<CoeffType>::TransportMixture( const ChemicalMixture<CoeffType>& chem_mix,
                                                 const std::string & filename, bool verbose, ParsingType type )
    : _chemical_mixture( chem_mix),
      _transport_species(_chemical_mixture.n_species(), NULL )
  {

    ParserBase<CoeffType> * parser(NULL);
    switch(type)
      {
      case ASCII:
        parser = new ASCIIParser<CoeffType>(filename,verbose);
        break;
      case CHEMKIN:
        parser = new ChemKinParser<CoeffType>(filename,verbose);
        break;
      case XML:
        parser = new XMLParser<CoeffType>(filename,verbose);
        break;
      default:
        antioch_parsing_error("unknown type");
      }

    // Now read in transport properties for the requested species and stash
    read_transport_species_data(parser,*this);

    // check we have everyone requested
    for( unsigned int s = 0; s < _transport_species.size(); ++s )
      {
        if(!_transport_species[s]) // it is not mandatory, Blottner or Sutherland are self-sufficient
          {
            std::cerr << "Warning: missing transport data for species "
                      << _chemical_mixture.species_inverse_name_map().at(_chemical_mixture.species_list()[s])
                      << "\n"
                      << "Be sure to use a transport model that does not require the default data"
                      << std::endl;
          }
      }

    delete parser;
  }

  template<typename CoeffType>
  TransportMixture<CoeffType>::TransportMixture( const ChemicalMixture<CoeffType> &mixture,
                                                 ParserBase<CoeffType> * parser)
    : _chemical_mixture( mixture),
      _transport_species(_chemical_mixture.n_species(), NULL )
  {
    // Now read in transport properties for the requested species and stash
    read_transport_species_data(parser,*this);

    // check we have everyone requested
    for( unsigned int s = 0; s < _transport_species.size(); ++s )
      {
        if(!_transport_species[s]) // it is not mandatory, Blottner or Sutherland are self-sufficient
          {
            std::cerr << "Warning: missing transport data for species " << _chemical_mixture.species_inverse_name_map().at(
                                                                                                                           _chemical_mixture.species_list()[s]) << "\n"
                      << "Be sure to use a transport model that does not require the default data"
                      << std::endl;
          }
      }

    return;

  }

  template<typename CoeffType>
  TransportMixture<CoeffType>::~TransportMixture()
  {
    // Clean up all the TransportSpecies we stored
    for( typename std::vector<TransportSpecies<CoeffType>* >::iterator it = _transport_species.begin();
         it < _transport_species.end(); ++it )
      {
        delete (*it);
      }
  }

  template<typename CoeffType>
  void TransportMixture<CoeffType>::add_species( const unsigned int index,
                                                 CoeffType LJ_depth,
                                                 CoeffType LJ_diameter,
                                                 CoeffType dipole_moment,
                                                 CoeffType polarizability,
                                                 CoeffType rotational_relaxation,
                                                 CoeffType mass)
  {
    Species name_enum = _chemical_mixture.species_list()[index];

    _transport_species[index] = new TransportSpecies<CoeffType>(name_enum,
                                                                LJ_depth,
                                                                LJ_diameter,
                                                                dipole_moment,
                                                                polarizability,
                                                                rotational_relaxation,mass);
  }

} // end namespace Antioch

// Instantiate
ANTIOCH_NUMERIC_TYPE_CLASS_INSTANTIATE(Antioch::TransportMixture);
