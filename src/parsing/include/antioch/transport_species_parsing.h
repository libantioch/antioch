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

#ifndef ANTIOCH_TRANSPORT_PARSING_H
#define ANTIOCH_TRANSPORT_PARSING_H

// Antioch
#include "antioch/ascii_parser.h"
#include "antioch/chemkin_parser.h"
#include "antioch/xml_parser.h"

//C++
#include <iostream>
#include <string>
#include <fstream>


namespace Antioch{

  //Forward declaration
  template <typename ThermoEvaluator, typename NumericType>
  class TransportMixture;


  template <typename ThermoEvaluator, typename NumericType>
  void read_transport_species_data(ParserBase<NumericType> * parser,
                                   TransportMixture<ThermoEvaluator,NumericType> & transport);

/*----------- inline functions ----------------*/


  template <typename ThermoEvaluator, typename NumericType>
  void read_transport_species_data_ascii(ParserBase<NumericType> * parser, TransportMixture<ThermoEvaluator,NumericType> & transport)
  {

    switch(parser->enum_type())
    {
      case ASCII:
      {
        (static_cast<ASCIIParser<NumericType> *>(parser))->read_transport_data(transport);
        break;
      }case CHEMKIN:
      {
        (static_cast<ChemKinParser<NumericType> *>(parser))->read_transport_data(transport);
        break;
      }case XML:
      {
        (static_cast<XMLParser<NumericType> *>(parser))->read_transport_data(transport);
        break;
      }default:
      {
        antioch_parsing_error("unknown parser type \"" + parser->type() + "\"!!!");
      }
    }

    // sanity check, we may require these informations
    bool fail(false);
    for(unsigned int s = 0; s < transport.chemical_species().size(); s++)
    {
        if(!transport.transport_species()[s])
        {
            fail = true;
            break;
        }
    }
    if(fail)
    {
      std::cerr << "Molecule(s) is(are) missing in transport description.  Please update the information."
                << "  Currently using file " << parser->file() << ".\n"
                << "You might have some problem later if you need these description.  "
                << "Missing molecule(s) is(are):" << std::endl;
      for(unsigned int i = 0; i < chem_mixture.species_list().size(); i++)
      {
        if(!transport.transport_species()[i])
        {
           std::cerr << transport.chem_mixture.species_inverse_name_map().at(i) << std::endl;
        }
      }
    }

    return;
  }

} //end namespace Antioch

#endif

