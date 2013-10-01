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

#ifndef ANTIOCH_REACTION_SET_DATA_XML_H
#define ANTIOCH_REACTION_SET_DATA_XML_H

// Antioch
#include "antioch/string_utils.h"
#include "antioch/reaction_set.h"
#include "antioch/kinetics_parsing.h"
#include "antioch/reaction_parsing.h"
#include "antioch/physical_constants.h"

// XML
#include "antioch/tinyxml2.h"

// C++
#include <string>
#include <vector>
#include <map>

namespace Antioch
{
  template<class NumericType>
  void read_reaction_set_data_xml( const std::string& filename,
                                   const bool verbose,
                                   ReactionSet<NumericType>& reaction_set );


  /* ------------------------- Inline Functions -------------------------*/
  template<class NumericType>
  inline
  void read_reaction_set_data_xml( const std::string& filename,
                                   const bool verbose,
                                   ReactionSet<NumericType>& reaction_set )
  {
    tinyxml2::XMLDocument doc;
    doc.LoadFile(filename.c_str());

    tinyxml2::XMLElement* element = doc.FirstChildElement("ctml");
    if (!element) 
      {
        std::cerr << "ERROR:  no <ctml> tag found in file " << filename
                  << std::endl;
        antioch_error();
      }
    
    const ChemicalMixture<NumericType>& chem_mixture = reaction_set.chemical_mixture();
    unsigned int n_species = chem_mixture.n_species();
    // Sanity Check on species
    /*
      tinyxml2::XMLElement* species = element->FirstChildElement("phase");
      species = species->FirstChildElement("speciesArray");

      std::vector<std::string> species_names; 
    
      std::cout << species->GetText() << std::endl;
    
      SplitString(std::string(species->GetText()),
      " ",
      species_names,
      false);

    
      if( n_species != species_names.size() )
      {
      std::cerr << "Error: Mismatch in n_species and the number of species specified in" << std::endl
      << "       the kinetics XML file " << filename << std::endl
      << "       Found species: " << species->GetText() << std::endl;
      antioch_error();
      }
    */

    // Now read in reaction data
    element = element->FirstChildElement("reactionData");
    if (!element) return;
    
    tinyxml2::XMLElement* reaction = element->FirstChildElement("reaction");

    std::map<KineticsModel::KineticsModel,std::string> kin_keyword;
    kin_keyword[KineticsModel::HERCOURT_ESSEN] = "HercourtEssen";
    kin_keyword[KineticsModel::BERTHELOT]      = "Berthelot";
    kin_keyword[KineticsModel::ARRHENIUS]      = "Arrhenius";
    kin_keyword[KineticsModel::BHE]            = "BerthelotHercourtEssen";
    kin_keyword[KineticsModel::KOOIJ]          = "Kooij";
    kin_keyword[KineticsModel::VANTHOFF]       = "VantHoff";

    std::vector<KineticsModel::KineticsModel> models; 
    models.push_back(KineticsModel::HERCOURT_ESSEN);
    models.push_back(KineticsModel::BERTHELOT);
    models.push_back(KineticsModel::ARRHENIUS);
    models.push_back(KineticsModel::BHE);
    models.push_back(KineticsModel::KOOIJ);
    models.push_back(KineticsModel::VANTHOFF);

    while (reaction)
      {
        if (verbose) std::cout << "Reaction #" << reaction->IntAttribute("id") << ":\n"
                               << " eqn: " << reaction->FirstChildElement("equation")->GetText()
                               << std::endl;

        ReactionType::ReactionType typeReaction(ReactionType::ELEMENTARY);
        KineticsModel::KineticsModel kineticsModel(KineticsModel::HERCOURT_ESSEN); // = 0

        if (reaction->Attribute("type"))
          {
            if (verbose) std::cout << " type: " << reaction->Attribute("type");
            if (std::string(reaction->Attribute("type")) == "threeBody")
              typeReaction = ReactionType::THREE_BODY;
          }
            
        // construct a Reaction object    
        Reaction<NumericType>* my_rxn = build_reaction<NumericType>(n_species, reaction->FirstChildElement("equation")->GetText(),typeReaction,kineticsModel);

        tinyxml2::XMLElement* rate_constant = reaction->FirstChildElement("rateCoeff")->FirstChildElement(kin_keyword[kineticsModel].c_str());
        unsigned int imod(0);
        while(!rate_constant)
        {
          if(imod == models.size())antioch_not_implemented();
          imod++;
          rate_constant = reaction->FirstChildElement("rateCoeff")->FirstChildElement(kin_keyword[models[imod]].c_str());
        }
        kineticsModel = models[imod];

        // usually Kooij is called Arrhenius, check here
        if(kineticsModel == KineticsModel::ARRHENIUS && rate_constant->FirstChildElement("b") != NULL)
        {
          if(std::atof(rate_constant->FirstChildElement("b")->GetText()) != 0.)
          {
              kineticsModel = KineticsModel::KOOIJ;
          }
        }

        if(verbose) 
          {
            std::cout << "\n rates: " << kin_keyword[kineticsModel] << " model\n"
                      << "   A: " << rate_constant->FirstChildElement("A")->GetText() << "\n"; //always
            if(rate_constant->FirstChildElement("b") != NULL)std::cout << "   b: " << rate_constant->FirstChildElement("b")->GetText() << "\n";
            if(rate_constant->FirstChildElement("E") != NULL)std::cout << "   E: " << rate_constant->FirstChildElement("E")->GetText() << "\n";
            if(rate_constant->FirstChildElement("D") != NULL)std::cout << "   D: " << rate_constant->FirstChildElement("D")->GetText() << "\n";
          }

        // typically Cantera files list activation energy in cal/mol, but we want it in K.
        std::vector<NumericType> data;
        data.push_back(std::atof(rate_constant->FirstChildElement("A")->GetText()));
        if(rate_constant->FirstChildElement("b") != NULL)data.push_back(std::atof(rate_constant->FirstChildElement("b")->GetText()));
        if(data.back() == 0.)data.pop_back();//if ARRHENIUS parameterized as KOOIJ
        if(rate_constant->FirstChildElement("E") != NULL)data.push_back(std::atof(rate_constant->FirstChildElement("E")->GetText()));
        if(rate_constant->FirstChildElement("D") != NULL)data.push_back(std::atof(rate_constant->FirstChildElement("D")->GetText()));
        //Tref
        if(kineticsModel == KineticsModel::HERCOURT_ESSEN ||
           kineticsModel == KineticsModel::BHE            ||
           kineticsModel == KineticsModel::KOOIJ          ||
           kineticsModel == KineticsModel::VANTHOFF) 
        {
          data.push_back(1.);
          if(rate_constant->FirstChildElement("Tref"))data.back() = std::atof(rate_constant->FirstChildElement("Tref")->GetText());
        }
        //scale E -> E/R
        if(kineticsModel == KineticsModel::ARRHENIUS ||
           kineticsModel == KineticsModel::KOOIJ     ||
           kineticsModel == KineticsModel::VANTHOFF)
        {
          data.push_back(Constants::R_universal<NumericType>()/1000.L);
          if( std::string(rate_constant->FirstChildElement("E")->Attribute("units")) == "cal/mol" )
            {
              data.back() = 1.9858775L;
            }
        }

        KineticsType<NumericType>* rate = build_rate<NumericType>(data,kineticsModel);

        my_rxn->add_forward_rate(rate);

        tinyxml2::XMLElement *efficiencies = 
          reaction->FirstChildElement("rateCoeff")->FirstChildElement("efficiencies");

        if(efficiencies)
          {
            if(efficiencies->GetText())
              {
                antioch_assert_equal_to (ReactionType::THREE_BODY, my_rxn->type());

                if(verbose) std::cout << "   efficiencies: " << efficiencies->GetText();

                std::vector<std::string> efficiency_pairs;

                SplitString( std::string(efficiencies->GetText()),
                             " ",
                             efficiency_pairs,
                             /* include_empties = */ false );

                for(unsigned int p = 0; p < efficiency_pairs.size(); p++)
                  {
                    std::pair<std::string, double> pair(split_string_double_on_colon (efficiency_pairs[p]));
                    if(verbose) 
                      {
                        std::cout  << "\n    " << efficiency_pairs[p] 
                                   << " " << pair.first << " " << pair.second;
                      }

                    if(pair.first == "e-") pair.first = "e";

                    // it is possible that the efficiency is specified for a species we are not
                    // modeling - so only add the efficiency if it is included in our list
                    if( chem_mixture.active_species_name_map().count( pair.first ) )
                      {
                        my_rxn->set_efficiency( pair.first,
                                               chem_mixture.active_species_name_map().find( pair.first )->second,
                                               pair.second );
                      }
                  }
              }
          }
        

        tinyxml2::XMLElement* reactants = reaction->FirstChildElement("reactants");
        tinyxml2::XMLElement* products  = reaction->FirstChildElement("products");
        
        // We will add the reaction, unless we do not have a 
        // reactant or product
        bool relevant_reaction = true;

        if(reactants->GetText())
          {
            if (verbose) std::cout << "\n   reactants: " << reactants->GetText();
                
            std::vector<std::string> reactant_pairs;
                
            // Split the reactant string on whitespace. If no entries were found,
            // there is no whitespace - and assume then only one reactant is listed.
            if( !SplitString(std::string(reactants->GetText()),
                             " ",
                             reactant_pairs,
                             /* include_empties = */ false)     )
              {
                reactant_pairs.push_back(reactants->GetText());
              }

            for( unsigned int p=0; p < reactant_pairs.size(); p++ )
              {
                std::pair<std::string,int> pair( split_string_int_on_colon(reactant_pairs[p]) );

                if(pair.first == "e-") pair.first = "e";

                if(verbose) std::cout  << "\n    " << reactant_pairs[p] << " " << pair.first << " " << pair.second;

                if( !chem_mixture.active_species_name_map().count( pair.first ) )
                  {
                    relevant_reaction = false;
                    if (verbose) std::cout << "\n     -> skipping this reaction (no reactant " << pair.first << ")";
                  }
                else
                  {
                    my_rxn->add_reactant( pair.first,
                                         chem_mixture.active_species_name_map().find( pair.first )->second,
                                         pair.second );
                  }
              }
          }
        if(products->GetText())
          {
            if(verbose) std::cout << "\n   products: " << products->GetText();
                
            std::vector<std::string> product_pairs;
                
            // Split the product string on whitespace. If no entries were found,
            // there is no whitespace - and assume then only one product is listed.
            if( !SplitString( std::string(products->GetText()),
                              " ",
                              product_pairs,
                              /* include_empties = */ false )     )
              {
                product_pairs.push_back(products->GetText());
              }

            for (unsigned int p=0; p<product_pairs.size(); p++)
              {
                std::pair<std::string, int> pair(split_string_int_on_colon (product_pairs[p]));

                if(pair.first == "e-") pair.first = "e";

                if(verbose) std::cout  << "\n    " << product_pairs[p] << " " << pair.first << " " << pair.second;

                if( !chem_mixture.active_species_name_map().count( pair.first ) )
                  {
                    relevant_reaction = false;
                    if (verbose) std::cout << "\n     -> skipping this reaction (no product " << pair.first << ")";
                  }
                else
                  {
                    my_rxn->add_product( pair.first,
                                        chem_mixture.active_species_name_map().find( pair.first )->second,
                                        pair.second );
                  }
              }
          }
        if(verbose) std::cout << "\n\n";

        if(relevant_reaction) 
          reaction_set.add_reaction(my_rxn);     

        // Go to the next reaction
        reaction = reaction->NextSiblingElement("reaction");
      }

    return;
  }
 
} // end namespace Antioch

#endif // ANTIOCH_REACTION_SET_DATA_XML_H
