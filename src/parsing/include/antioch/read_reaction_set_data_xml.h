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
    
    while (reaction)
      {
	if (verbose) std::cout << "Reaction #" << reaction->IntAttribute("id") << ":\n"
			       << " eqn: " << reaction->FirstChildElement("equation")->GetText()
			       << std::endl;

        ReactionType::ReactionType typeReaction(ReactionType::ELEMENTARY);
        KineticsModel::KineticsModel kineticsModel(KineticsModel::KOOIJ);

	if (reaction->Attribute("type"))
	  {
	    if (verbose) std::cout << " type: " << reaction->Attribute("type");
	    if (std::string(reaction->Attribute("type")) == "threeBody")
	      typeReaction = ReactionType::THREE_BODY;
	  }
	    
	// construct a Reaction object	  
	Reaction<NumericType>* my_rxn = build_reaction<NumericType>(n_species, reaction->FirstChildElement("equation")->GetText(),typeReaction,kineticsModel);

	tinyxml2::XMLElement* Arrhenius = reaction->FirstChildElement("rateCoeff")->FirstChildElement("Arrhenius");

	if(verbose) 
	  {
	    std::cout << "\n rates: Kooij model\n"
		      << "   A: " << Arrhenius->FirstChildElement("A")->GetText() << "\n"
		      << "   b: " << Arrhenius->FirstChildElement("b")->GetText() << "\n"
		      << "   E: " << Arrhenius->FirstChildElement("E")->GetText() << "\n";
	  }

	// typically Cantera files list activation energy in cal/mol, but we want it in K.
        std::vector<NumericType> data;
        data.push_back(std::atof(Arrhenius->FirstChildElement("A")->GetText()));
        data.push_back(std::atof(Arrhenius->FirstChildElement("b")->GetText()));
        data.push_back(std::atof(Arrhenius->FirstChildElement("E")->GetText()));
        data.push_back(1.);
        data.push_back(Constants::R_universal<NumericType>()/1000.);
        
        if(Arrhenius->FirstChildElement("Tref"))data[3] = std::atof(Arrhenius->FirstChildElement("Tref")->GetText());
	if( std::string(Arrhenius->FirstChildElement("E")->Attribute("units")) == "cal/mol" )
	  {
	    data[4] = 1.9858775;
	  }
	else 
	  {
	    // huh?
	    antioch_error();
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
