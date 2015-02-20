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

#ifndef ANTIOCH_SPECIES_ASCII_PARSING_H
#define ANTIOCH_SPECIES_ASCII_PARSING_H

// Antioch
#include "antioch/ascii_parser.h"
#include "antioch/chemkin_parser.h"
#include "antioch/xml_parser.h"

// C++
#include <fstream>

namespace Antioch
{
  typedef unsigned int Species;

  // Forward declarations
  template <class NumericType>
  class ChemicalMixture;

    
  template <typename NumericType>
  void read_chemical_species_composition(ParserBase<NumericType> * parser,
                                         ChemicalMixture<NumericType> & mixture);

  template<class NumericType>
  void read_species_data( ParserBase<NumericType> * parser,
                          ChemicalMixture<NumericType>& chem_mixture);
				

  template<class NumericType>
  void read_species_vibrational_data(ParserBase<NumericType> * parser,
                                     ChemicalMixture<NumericType>& chem_mixture);

  template<class NumericType>
  void read_species_electronic_data(ParserBase<NumericType> * parser,
                                    ChemicalMixture<NumericType>& chem_mixture);

//----------------------------------------------------------------

   template <typename NumericType>
   inline
   void read_chemical_species_composition(ParserBase<NumericType> * parser,
                                          ChemicalMixture<NumericType> & mixture)
   {
      std::vector<std::string> species;
      switch(parser->enum_type())
      {
        case ASCII:
        {
           species = (static_cast<ASCIIParser<NumericType> *>(parser))->species_list();
           break;
        }case CHEMKIN:
        {
           species = static_cast<ChemKinParser<NumericType> *>(parser)->species_list();
           break;
        }case XML:
        {
           species = static_cast<XMLParser<NumericType> *>(parser)->species_list();
           break;
        }default:
        {
           antioch_parsing_error("unknown parser type \"" + parser->type() + "\"!!!");
        }
      }
      mixture.initialize_species(species);
   }
  
  template<class NumericType>
  inline
  void read_species_data(ParserBase<NumericType> * parser,
                         ChemicalMixture<NumericType>& chem_mixture)
  {
    switch(parser->enum_type())
    {
      case ASCII:
      {
        (static_cast<ASCIIParser<NumericType> *>(parser))->read_chemical_species(chem_mixture);
        break;
      }case CHEMKIN:
      {
        (static_cast<ChemKinParser<NumericType> *>(parser))->read_chemical_species(chem_mixture);
        break;
      }case XML:
      {
        (static_cast<XMLParser<NumericType> *>(parser))->read_chemical_species(chem_mixture);
        break;
      }default:
      {
        antioch_parsing_error("unknown parser type \"" + parser->type() + "\"!!!");
      }
    }

    // sanity check, we require these informations
    bool fail(false);
    for(unsigned int s = 0; s < chem_mixture.chemical_species().size(); s++)
    {
        if(!chem_mixture.chemical_species()[s])
        {
            fail = true;
            break;
        }
    }
    if(fail)
    {
      std::cerr << "Molecule(s) is(are) missing.  Please update the information."
                << "  Currently using file " << parser->file() << ".\n"
                << "Missing molecule(s) is(are):" << std::endl;
      for(unsigned int i = 0; i < chem_mixture.species_list().size(); i++)
      {
        if(!chem_mixture.chemical_species()[i])
        {
           std::cerr << chem_mixture.species_inverse_name_map().at(i) << std::endl;
        }
      }
      antioch_error();
    }

    return;
  }

  template<class NumericType>
  inline
  void read_species_vibrational_data(ParserBase<NumericType> * parser,
                                     ChemicalMixture<NumericType>& chem_mixture)
  {
    switch(parser->enum_type())
    {
      case ASCII:
      {
        (static_cast<ASCIIParser<NumericType> *>(parser))->read_vibrational_data(chem_mixture);
        break;
      }case CHEMKIN:
      {
        (static_cast<ChemKinParser<NumericType> *>(parser))->read_vibrational_data(chem_mixture);
        break;
      }case XML:
      {
        (static_cast<XMLParser<NumericType> *>(parser))->read_vibrational_data(chem_mixture);
        break;
      }default:
      {
        antioch_parsing_error("unknown parser type \"" + parser->type() + "\"!!!");
      }
    }


    // sanity check, we check these informations
    std::vector<std::string> missing;
    for(unsigned int s = 0; s < chem_mixture.chemical_species().size(); s++)
    {
        if(chem_mixture.chemical_species()[s]->theta_v().empty())missing.push_back(chem_mixture.chemical_species()[s]->species());
    }
    if(!missing.empty())
    {
       std::cerr << "WARNING:\nVibrational levels are missing.  Please update the information."
                 << "  Currently using file " << parser->file() << ".\n"
                 << "Missing molecule(s) is(are):" << std::endl;
       for(unsigned int m = 0; m < missing.size(); m++)std::cerr << missing[m] << std::endl;
    }
    return;
  }

  
  template<class NumericType>
  inline
  void read_species_electronic_data(ParserBase<NumericType> * parser,
                                    ChemicalMixture<NumericType>& chem_mixture)
                                           
  {
    switch(parser->enum_type())
    {
      case ASCII:
      {
        (static_cast<ASCIIParser<NumericType> *>(parser))->read_electronic_data(chem_mixture);
        break;
      }case CHEMKIN:
      {
        (static_cast<ChemKinParser<NumericType> *>(parser))->read_electronic_data(chem_mixture);
        break;
      }case XML:
      {
        (static_cast<XMLParser<NumericType> *>(parser))->read_electronic_data(chem_mixture);
        break;
      }default:
      {
        antioch_parsing_error("unknown parser type \"" + parser->type() + "\"!!!");
      }
    }
    
    // sanity check, we check these informations
    std::vector<std::string> missing;
    for(unsigned int s = 0; s < chem_mixture.chemical_species().size(); s++)
    {
        if(chem_mixture.chemical_species()[s]->theta_e().empty())missing.push_back(chem_mixture.chemical_species()[s]->species());
    }
    if(!missing.empty())
    {
       std::cerr << "WARNING:\nElectronic levels are missing.  Please update the information."
                 << "  Currently using file " << parser->file() << ".\n"
                 << "Missing molecule(s) is(are):" << std::endl;
       for(unsigned int m = 0; m < missing.size(); m++)std::cerr << missing[m] << std::endl;
    }
    return;
  }
  
} // end namespace Antioch

#endif // ANTIOCH_SPECIES_ASCII_PARSING_H
