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

#include "antioch/species_parsing.h"

// Antioch
#include "antioch/parser_base.h"
#include "antioch/chemical_mixture.h"
#include "antioch/species_parsing_instantiation_macro.h"

// C++
#include <iostream>

namespace Antioch
{
  template <typename NumericType>
  void read_chemical_species_composition(ParserBase<NumericType> * parser,
                                         ChemicalMixture<NumericType> & mixture)
  {
    std::vector<std::string> species = parser->species_list();

    mixture.initialize_species(species);
  }

  template<class NumericType>
  void read_species_data(ParserBase<NumericType> * parser,
                         ChemicalMixture<NumericType>& chem_mixture)
  {
    parser->read_chemical_species(chem_mixture);

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
        std::stringstream ss;
        ss << "Molecule(s) is(are) missing.  Please update the information."
           << "  Currently using file " << parser->file() << ".\n"
           << "Missing molecule(s) is(are):" << std::endl;
        for(unsigned int i = 0; i < chem_mixture.species_list().size(); i++)
          {
            if(!chem_mixture.chemical_species()[i])
              {
                ss << chem_mixture.species_inverse_name_map().at(i) << std::endl;
              }
          }
        antioch_error_msg(ss.str());
      }
  }

  template<class NumericType>
  void read_species_vibrational_data(ParserBase<NumericType> * parser,
                                     ChemicalMixture<NumericType>& chem_mixture)
  {
    parser->read_vibrational_data(chem_mixture);

    // sanity check, we check these informations
    std::vector<std::string> missing;
    for(unsigned int s = 0; s < chem_mixture.chemical_species().size(); s++)
      {
        if(chem_mixture.chemical_species()[s]->theta_v().empty())missing.push_back(chem_mixture.chemical_species()[s]->species());
      }
    if(!missing.empty())
      {
        std::stringstream ss;
        ss << "WARNING:\nVibrational levels are missing.  Please update the information."
           << "  Currently using file " << parser->file() << ".\n"
           << "Missing molecule(s) is(are):" << std::endl;
        for(unsigned int m = 0; m < missing.size(); m++)
          ss << missing[m] << std::endl;

        antioch_warning(ss.str());
      }
  }


  template<class NumericType>
  void read_species_electronic_data(ParserBase<NumericType> * parser,
                                    ChemicalMixture<NumericType>& chem_mixture)

  {
    parser->read_electronic_data(chem_mixture);

    // sanity check, we check these informations
    std::vector<std::string> missing;
    for(unsigned int s = 0; s < chem_mixture.chemical_species().size(); s++)
      {
        if(chem_mixture.chemical_species()[s]->theta_e().empty())missing.push_back(chem_mixture.chemical_species()[s]->species());
      }
    if(!missing.empty())
      {
        std::stringstream ss;
        ss << "WARNING:\nElectronic levels are missing.  Please update the information."
           << "  Currently using file " << parser->file() << ".\n"
           << "Missing molecule(s) is(are):" << std::endl;
        for(unsigned int m = 0; m < missing.size(); m++)
          ss << missing[m] << std::endl;

        antioch_warning(ss.str());
      }
  }

  ANTIOCH_SPECIES_PARSING_INSTANTIATE();

} // end namespace Antioch
