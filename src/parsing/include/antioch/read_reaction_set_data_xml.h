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
#include "antioch/antioch_asserts.h"
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
    if(doc.LoadFile(filename.c_str()))
      {
        std::cerr << "ERROR: unable to load xml file " << filename << std::endl;
        std::cerr << "Error of tinyxml2 library:\n"
                  << "\tID = " << doc.ErrorID() << "\n"
                  << "\tError String1 = " << doc.GetErrorStr1() << "\n"
                  << "\tError String2 = " << doc.GetErrorStr2() << std::endl;
        antioch_error();
      }


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

    std::map<std::string,KineticsModel::KineticsModel> kin_keyword;
    kin_keyword["HercourtEssen"]          = KineticsModel::HERCOURT_ESSEN;
    kin_keyword["Berthelot"]              = KineticsModel::BERTHELOT;
    kin_keyword["Arrhenius"]              = KineticsModel::ARRHENIUS;
    kin_keyword["BerthelotHercourtEssen"] = KineticsModel::BHE;
    kin_keyword["Kooij"]                  = KineticsModel::KOOIJ;
    kin_keyword["ModifiedArrhenius"]      = KineticsModel::KOOIJ;  //for Arrhenius fans
    kin_keyword["VantHoff"]               = KineticsModel::VANTHOFF;

    std::vector<std::string> models; 
    models.push_back("HercourtEssen");
    models.push_back("Berthelot");
    models.push_back("Arrhenius");
    models.push_back("BerthelotHercourtEssen");
    models.push_back("Kooij");
    models.push_back("ModifiedArrhenius");
    models.push_back("VantHoff");

    std::map<std::string,ReactionType::ReactionType> proc_keyword;
    proc_keyword["Elementary"]        = ReactionType::ELEMENTARY;
    proc_keyword["Duplicate"]         = ReactionType::DUPLICATE;
    proc_keyword["ThreeBody"]         = ReactionType::THREE_BODY;
    proc_keyword["LindemannFalloff"]  = ReactionType::LINDEMANN_FALLOFF;
    proc_keyword["TroeFalloff"]       = ReactionType::TROE_FALLOFF;

    while (reaction)
      {
        if (verbose) std::cout << "Reaction #" << reaction->IntAttribute("id") << ":\n"
                               << " eqn: " << reaction->FirstChildElement("equation")->GetText()
                               << std::endl;

        ReactionType::ReactionType typeReaction(ReactionType::ELEMENTARY);
        KineticsModel::KineticsModel kineticsModel(KineticsModel::HERCOURT_ESSEN); // = 0

        if (reaction->Attribute("type"))
          {
            if (verbose) std::cout << " type: " << reaction->Attribute("type") << std::endl;
            if(!proc_keyword.count(reaction->Attribute("type")))
            {
                std::cerr << "Implemented chemical processes are:\n"
                          << "  Elementary (default)\n"
                          << "  Duplicate\n"
                          << "  ThreeBody\n"
                          << "  LindemannFalloff\n"
                          << "  TroeFalloff\n" 
                          << "See Antioch documentation for more details."
                          << std::endl;
                antioch_not_implemented();
            }
            typeReaction = proc_keyword[reaction->Attribute("type")];
          }
            
        unsigned int imod(0);
        tinyxml2::XMLElement* rate_constant = reaction->FirstChildElement("rateCoeff")->FirstChildElement(models[imod].c_str());
        while(!rate_constant)
        {
          if(imod == models.size() - 1)
          {
                std::cerr << "Could not find a suitable kinetics model.\n"
                          << "Implemented kinetics models are:\n"
                          << "  HercourtEssen\n"
                          << "  Berthelot\n"
                          << "  Arrhenius\n"
                          << "  BerthelotHercourtEssen\n"
                          << "  Kooij (default)\n"
                          << "  ModifiedArrhenius (= Kooij)\n"
                          << "  VantHoff\n"
                          << "See Antioch documentation for more details."
                          << std::endl;
                antioch_not_implemented();
          }
          imod++;
          rate_constant = reaction->FirstChildElement("rateCoeff")->FirstChildElement(models[imod].c_str());
        }
        kineticsModel = kin_keyword[models[imod]];

        // construct a Reaction object    
        Reaction<NumericType>* my_rxn = build_reaction<NumericType>(n_species, reaction->FirstChildElement("equation")->GetText(),typeReaction,kineticsModel);

        while(rate_constant) //for duplicate and falloff models, several kinetics rate to load, no mixing allowed
        {

          // usually Kooij is called Arrhenius, check here
          if(kineticsModel == KineticsModel::ARRHENIUS && rate_constant->FirstChildElement("b") != NULL)
          {
            if(std::atof(rate_constant->FirstChildElement("b")->GetText()) != 0.)
            {
               kineticsModel = KineticsModel::KOOIJ;
               antioch_deprecated();
               std::cerr << "An equation of the form \"A * (T/Tref)^beta * exp(-Ea/(R*T))\" is a Kooij equation,\n"
                         << "I guess a modified Arrhenius could be a name too.  Whatever, the correct label is\n"
                         << "\"Kooij\", or, << à la limite >> \"ModifiedArrhenius\".  Please use those terms instead,\n"
                         << "thanks and a good day to you, user." << std::endl;
            }
          }

          if(verbose) 
            {
              std::cout << " rate: " << models[imod] << " model\n"
                        << "   A: " << rate_constant->FirstChildElement("A")->GetText() << "\n"; //always
              if(rate_constant->FirstChildElement("b") != NULL)
              {
                  std::cout << "   b: " << rate_constant->FirstChildElement("b")->GetText() << "\n";
              }
              if(rate_constant->FirstChildElement("E") != NULL)
              {
                  std::cout << "   E: " << rate_constant->FirstChildElement("E")->GetText() << "\n";
              }
              if(rate_constant->FirstChildElement("D") != NULL)
              {
                  std::cout << "   D: " << rate_constant->FirstChildElement("D")->GetText() << "\n";
              }
            }

          // typically Cantera files list activation energy in cal/mol, but we want it in K.
          std::vector<NumericType> data;
          data.push_back(std::atof(rate_constant->FirstChildElement("A")->GetText()));
          if(rate_constant->FirstChildElement("b") != NULL)
          {
             data.push_back(std::atof(rate_constant->FirstChildElement("b")->GetText()));
          }
          if(data.back() == 0.)//if ARRHENIUS parameterized as KOOIJ
          {
             data.pop_back();
          }
          if(rate_constant->FirstChildElement("E") != NULL)
          {
             data.push_back(std::atof(rate_constant->FirstChildElement("E")->GetText()));
          }
          if(rate_constant->FirstChildElement("D") != NULL)
          {
             data.push_back(std::atof(rate_constant->FirstChildElement("D")->GetText()));
          }
          //Tref
          if(kineticsModel == KineticsModel::HERCOURT_ESSEN ||
             kineticsModel == KineticsModel::BHE            ||
             kineticsModel == KineticsModel::KOOIJ          ||
             kineticsModel == KineticsModel::VANTHOFF) 
          {
            data.push_back(1.);
            if(rate_constant->FirstChildElement("Tref"))
            {
                data.back() = std::atof(rate_constant->FirstChildElement("Tref")->GetText());
            }
          }
          //scale E -> E/R
          if(kineticsModel == KineticsModel::ARRHENIUS ||
             kineticsModel == KineticsModel::KOOIJ     ||
             kineticsModel == KineticsModel::VANTHOFF)
          {
            data.push_back(Constants::R_universal<NumericType>()/1000.L);
            if( rate_constant->FirstChildElement("E")->Attribute("units"))//if there's the attribute
              {
              if( std::string(rate_constant->FirstChildElement("E")->Attribute("units")) == "cal/mol" )
                {
                  data.back() = 1.9858775L;
                }
              }
          }

          KineticsType<NumericType>* rate = build_rate<NumericType>(data,kineticsModel);

          my_rxn->add_forward_rate(rate);

          rate_constant = rate_constant->NextSiblingElement(models[imod].c_str());

        } //end of duplicate/falloff kinetics description loop

        // for falloff, we need a way to know which rate constant is the low pressure limit 
        // and which is the high pressure limit
        // usually by calling the low pressure limite "k0". If nothing given, by default
        // the first rate constant encountered is the low limit,
        // so we need to change something only if the second rate constant has a "name" attribute
        // of value "k0"
        if(typeReaction == ReactionType::LINDEMANN_FALLOFF ||
           typeReaction == ReactionType::TROE_FALLOFF)
        {
           antioch_assert_equal_to(my_rxn->n_rate_constants(),2);
           rate_constant = reaction->FirstChildElement("rateCoeff")->FirstChildElement(models[imod].c_str())->NextSiblingElement(models[imod].c_str());
           if(rate_constant->Attribute("name")) //if attribute exists
           { 
             if(std::string(rate_constant->Attribute("name")) == "k0") //and is indeed k0
             {
               my_rxn->swap_forward_rates(0,1);
             }
           }
        }

        tinyxml2::XMLElement *efficiencies = 
          reaction->FirstChildElement("rateCoeff")->FirstChildElement("efficiencies");

        //efficiencies are only for three body reactions
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

        //F parameters only for Troe falloff
        tinyxml2::XMLElement *troe = 
          reaction->FirstChildElement("rateCoeff")->FirstChildElement("Troe");
        if(troe)
        {
std::cout << "0" << std::endl;
           antioch_assert_equal_to (ReactionType::TROE_FALLOFF, my_rxn->type());
           FalloffReaction<NumericType,TroeFalloff<NumericType> > *my_fall_rxn =
                static_cast<FalloffReaction<NumericType,TroeFalloff<NumericType> > *> (my_rxn);

           if(!troe->FirstChildElement("alpha"))
           {
                std::cerr << "alpha parameter of Troe falloff missing!" << std::endl;
                antioch_error();
           }
           my_fall_rxn->F().set_alpha(std::atof(troe->FirstChildElement("alpha")->GetText()));
           if(!troe->FirstChildElement("T3"))
           {
                std::cerr << "T*** parameter of Troe falloff missing!" << std::endl;
                antioch_error();
           }
           my_fall_rxn->F().set_T3(std::atof(troe->FirstChildElement("T3")->GetText()));
           if(!troe->FirstChildElement("T1"))
           {
                std::cerr << "T* parameter of Troe falloff missing!" << std::endl;
                antioch_error();
           }
           my_fall_rxn->F().set_T1(std::atof(troe->FirstChildElement("T1")->GetText()));
           if(troe->FirstChildElement("T2"))//T2 is optional
           {
             my_fall_rxn->F().set_T2(std::atof(troe->FirstChildElement("T2")->GetText()));
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
