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

#ifndef ANTIOCH_REACTION_SET_DATA_XML_H
#define ANTIOCH_REACTION_SET_DATA_XML_H

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/string_utils.h"
#include "antioch/reaction_set.h"
#include "antioch/kinetics_parsing.h"
#include "antioch/reaction_parsing.h"
#include "antioch/physical_constants.h"
#include "antioch/units.h"

// XML
#include "antioch/tinyxml2_imp.h"

// C++
#include <string>
#include <vector>
#include <map>

namespace Antioch
{

 /*!\file read_reaction_set_data_xml.h
  *
  * We parse the XML file here, with an exhaustive
  * unit management. The starting point is the kinetics
  * equation:
  * \f[
  *     \frac{\partial c}{\partial t} = k \prod_{s \in \text{reactants}} c_s^{l_s}
  * \f]
  * with \f$l\f$ the partial order of the reaction with respect to reactants \f$s\f$.
  * We obtain thus
  * \f[
  *     \unit{[k] = [M]^{m - 1} s^{-1}}
  * \f]
  * with \f$m\f$ the order of the reaction. By definition
  * \f[
  *     m = \sum_{s \in \text{reactants}} l_s
  * \f]
  * We are in an elementary processes paradigm, thus for a reactant species \f$s\f$,
  * \f$l_s = -\nu_s\f$ with \f$\nu_s\f$ the stoichiometric coefficient of reactant
  * species \f$s\f$.
  *
  * Example:
  * \f[
  *   \begin{array}{c}
  *      \ce{a A + b B -> c C + d D} \\
  *      m = a + b
  *   \end{array}
  * \f]
  *
  * To this, we consider the kinetics model (they're all included in the
  * Van't Hoff equation):
  * \f[
  *   \alpha(T) = A \left(\frac{T}{\mathrm{T_\text{ref}}}\right)^\beta\exp\left(-\frac{E_a}{\mathrm{R}T} + D T\right)
  * \f]
  *
  * We derive from this all the tests and default units:
  * \f[
  *  \begin{array}{lcccc}\toprule
  *                                & A                                           & \beta & E_a                & D \\\midrule
  *   \text{Elementary}            & \unit{\left(m^3mol^{-1}\right)^{m-1}s^{-1}} & -     & \unit{J\,mol^{-1}} & \unit{K^{-1}} \\
  *   \text{Duplicate}             & \unit{\left(m^3mol^{-1}\right)^{m-1}s^{-1}} & -     & \unit{J\,mol^{-1}} & \unit{K^{-1}} \\
  *   \text{Three body}            & \unit{\left(m^3mol^{-1}\right)^{m}s^{-1}}   & -     & \unit{J\,mol^{-1}} & \unit{K^{-1}} \\
  *   \text{Falloff}\quad k_0      & \unit{\left(m^3mol^{-1}\right)^{m}s^{-1}}   & -     & \unit{J\,mol^{-1}} & \unit{K^{-1}} \\
  *   \text{Falloff}\quad k_\infty & \unit{\left(m^3mol^{-1}\right)^{m-1}s^{-1}} & -     & \unit{J\,mol^{-1}} & \unit{K^{-1}} \\\bottomrule
  *  \end{array}
  * \f]
  * for the Troe falloff, the additionnal parameters are:
  * \f[
  *  \begin{array}{cccc}\toprule
  *    \alpha  & T^*      & T^{**}   & T^{***} \\\midrule
  *     -      & \unit{K} & \unit{K} & \unit{K} \\\bottomrule
  *  \end{array}
  * \f]
  *
  * Thus the reading is made in this fashion:
  *   - read reactants and products, get \f$m\f$
  *   - find default unit of \f$A\f$
  *   - read other parameters
  */
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
    kin_keyword["Constant"]               = KineticsModel::CONSTANT;
    kin_keyword["HercourtEssen"]          = KineticsModel::HERCOURT_ESSEN;
    kin_keyword["Berthelot"]              = KineticsModel::BERTHELOT;
    kin_keyword["Arrhenius"]              = KineticsModel::ARRHENIUS;
    kin_keyword["BerthelotHercourtEssen"] = KineticsModel::BHE;
    kin_keyword["Kooij"]                  = KineticsModel::KOOIJ;
    kin_keyword["ModifiedArrhenius"]      = KineticsModel::KOOIJ;  //for Arrhenius fans
    kin_keyword["VantHoff"]               = KineticsModel::VANTHOFF;
    kin_keyword["photochemistry"]         = KineticsModel::PHOTOCHEM;

    std::vector<std::string> models; 
    models.push_back("Constant");
    models.push_back("HercourtEssen");
    models.push_back("Berthelot");
    models.push_back("Arrhenius");
    models.push_back("BerthelotHercourtEssen");
    models.push_back("Kooij");
    models.push_back("ModifiedArrhenius");
    models.push_back("VantHoff");
    models.push_back("photochemistry");

    std::map<std::string,ReactionType::ReactionType> proc_keyword;
    proc_keyword["Elementary"]        = ReactionType::ELEMENTARY;
    proc_keyword["Duplicate"]         = ReactionType::DUPLICATE;
    proc_keyword["ThreeBody"]         = ReactionType::THREE_BODY;
    proc_keyword["threeBody"]         = ReactionType::THREE_BODY; // Cantera/backward compatiblity
    proc_keyword["LindemannFalloff"]  = ReactionType::LINDEMANN_FALLOFF;
    proc_keyword["TroeFalloff"]       = ReactionType::TROE_FALLOFF;

    while (reaction)
      {
        if (verbose) std::cout << "Reaction #" << reaction->IntAttribute("id") << ":\n"
                               << " eqn: " << reaction->FirstChildElement("equation")->GetText()
                               << std::endl;

        ReactionType::ReactionType typeReaction(ReactionType::ELEMENTARY);
        KineticsModel::KineticsModel kineticsModel(KineticsModel::HERCOURT_ESSEN); // = 0
        bool reversible(true);

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
            
        if(reaction->Attribute("reversible"))
        {
          if (verbose) std::cout << "reversible: " << reaction->Attribute("reversible") << std::endl;
          if(std::string(reaction->Attribute("reversible")) == "no")reversible = false;
        }

        unsigned int imod(0);
        tinyxml2::XMLElement* rate_constant = reaction->FirstChildElement("rateCoeff")->FirstChildElement(models[imod].c_str());
        while(!rate_constant)
        {
          if(imod == models.size() - 1)
          {
                std::cerr << "Could not find a suitable kinetics model.\n"
                          << "Implemented kinetics models are:\n"
                          << "  Constant\n"
                          << "  HercourtEssen\n"
                          << "  Berthelot\n"
                          << "  Arrhenius\n"
                          << "  BerthelotHercourtEssen\n"
                          << "  Kooij (default)\n"
                          << "  ModifiedArrhenius (= Kooij)\n"
                          << "  VantHoff\n"
                          << "  photochemistry\n"
                          << "See Antioch documentation for more details."
                          << std::endl;
                antioch_not_implemented();
          }
          imod++;
          rate_constant = reaction->FirstChildElement("rateCoeff")->FirstChildElement(models[imod].c_str());
        }
        kineticsModel = kin_keyword[models[imod]];

        tinyxml2::XMLElement* reactants = reaction->FirstChildElement("reactants");
        tinyxml2::XMLElement* products  = reaction->FirstChildElement("products");
        
        // usually Kooij is called Arrhenius, check here
        if(kineticsModel == KineticsModel::ARRHENIUS && rate_constant->FirstChildElement("b") != NULL)
        {
          if(rate_constant->FirstChildElement("b"))
          {
            if(std::atof(rate_constant->FirstChildElement("b")->GetText()) != 0.)
            {
               kineticsModel = KineticsModel::KOOIJ;
               std::cerr << "In reaction " << reaction->Attribute("id") << "\n"
                         << "An equation of the form \"A * (T/Tref)^beta * exp(-Ea/(R*T))\" is a Kooij equation,\n"
                         << "I guess a modified Arrhenius could be a name too.  Whatever, the correct label is\n"
                         << "\"Kooij\", or, << à la limite >> \"ModifiedArrhenius\".  Please use those terms instead,\n"
                         << "thanks and a good day to you, user." << std::endl;
            }
          }
        }


        // construct a Reaction object    
        Reaction<NumericType>* my_rxn = build_reaction<NumericType>(n_species, reaction->FirstChildElement("equation")->GetText(),
                                                                               reversible,typeReaction,kineticsModel);

        // We will add the reaction, unless we do not have a 
        // reactant or product
        bool relevant_reaction = true;
        unsigned int order_reaction(0);

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

                if( !chem_mixture.species_name_map().count( pair.first ) )
                  {
                    relevant_reaction = false;
                    if (verbose) std::cout << "\n     -> skipping this reaction (no reactant " << pair.first << ")";
                  }
                else
                  {
                    my_rxn->add_reactant( pair.first,
                                         chem_mixture.species_name_map().find( pair.first )->second,
                                         pair.second );
                    order_reaction += pair.second;
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

                if( !chem_mixture.species_name_map().count( pair.first ) )
                  {
                    relevant_reaction = false;
                    if (verbose) std::cout << "\n     -> skipping this reaction (no product " << pair.first << ")";
                  }
                else
                  {
                    my_rxn->add_product( pair.first,
                                        chem_mixture.species_name_map().find( pair.first )->second,
                                        pair.second );
                  }
              }
          }


        if(!relevant_reaction)
        {
          if(verbose) std::cout << "skipped reaction\n\n";
          reaction = reaction->NextSiblingElement("reaction");
          delete my_rxn;
          continue;
        }

        while(rate_constant) //for duplicate and falloff models, several kinetics rate to load, no mixing allowed
        {
        // typically Cantera files list 
        //      activation energy in cal/mol, but we want it in K.
        //      pre-exponential parameters in (m3/kmol)^(m-1)/s
        //      power parameter without unit
        // if falloff, we need to know who's k0 and kinfty
        // if photochemistry, we have a cross-section on a lambda grid
        //      cross-section typically in cm2/nm (cross-section on a resolution bin, 
        //                                          if bin unit not given, it is lambda unit (supposed to anyway), and a warning message)
        //      lambda typically in nm, sometimes in ang, default considered here is nm
        //                         you can also have cm-1, conversion is done with
        //                         formulae nm = cm-1 * / * adapted factor
          std::vector<NumericType> data;
          Units<NumericType> def_unit;
          int pow_unit(order_reaction - 1);
         //threebody always
          if(my_rxn->type() == ReactionType::THREE_BODY)pow_unit++;
        //falloff for k0
          if(my_rxn->type() == ReactionType::LINDEMANN_FALLOFF ||
             my_rxn->type() == ReactionType::TROE_FALLOFF)
          {
//k0 is either determined by an explicit name, or is the first of unnamed rate constants
             if(rate_constant->Attribute("name"))
             {
                if(std::string(rate_constant->Attribute("name")) == "k0")pow_unit++;
             }else if(my_rxn->n_rate_constants() == 0) // if we're indeed at the first reading
             {
                if(!rate_constant->NextSiblingElement(models[imod].c_str())->Attribute("name")) // and the next doesn't have a name
                     pow_unit++;
             }
          }
          def_unit.set_unit("m3/kmol"); //default
          def_unit *= pow_unit; //to the m-1 power
          def_unit.substract("s"); // per second
          
          if(verbose) 
            {
              std::cout << " rate: " << models[imod] << " model\n" << "\n";
              if(rate_constant->FirstChildElement("A") != NULL)
              {
                  std::cout << "   A: " << rate_constant->FirstChildElement("A")->GetText() << "\n";
              }
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
              if(rate_constant->FirstChildElement("lambda") != NULL)
              {
                  std::cout << "   lambda:\n" << rate_constant->FirstChildElement("lambda")->GetText() << "\n";
              }
              if(rate_constant->FirstChildElement("cross_section") != NULL)
              {
                  std::cout << "   cross section:\n" << rate_constant->FirstChildElement("cross_section")->GetText() << "\n";
              }
            }

          if(rate_constant->FirstChildElement("A") != NULL)
          {
            if(!rate_constant->FirstChildElement("A")->Attribute("units"))
            {
               antioch_unit_required("A",def_unit.get_symbol());
            }else
            {
              Units<NumericType> read_unit;
              read_unit.set_unit(rate_constant->FirstChildElement("A")->Attribute("units"));
              if(!read_unit.is_homogeneous(def_unit))
              {
                std::string errorstring("Error in reaction " + my_rxn->equation());
                errorstring += "\n A units should be homogeneous to " + def_unit.get_symbol() + 
                               " and you provided " + read_unit.get_symbol();
                antioch_unit_error(errorstring);
              }
              def_unit = read_unit;
           }
           if(verbose) 
             {
             std::cout  << "   A: " << rate_constant->FirstChildElement("A")->GetText()
                        << " " << def_unit.get_symbol() << std::endl; 
             }
           data.push_back(std::atof(rate_constant->FirstChildElement("A")->GetText()) * def_unit.get_SI_factor());
         }
//b has no unit
         if(rate_constant->FirstChildElement("b") != NULL)
         {
            data.push_back(std::atof(rate_constant->FirstChildElement("b")->GetText()));
            if(data.back() == 0.)//if ARRHENIUS parameterized as KOOIJ
            {
               data.pop_back();
               std::cerr << "In reaction " << reaction->Attribute("id") << "\n"
                         << "An equation of the form \"A * exp(-Ea/(R*T))\" is an Arrhenius equation,\n"
                         << "and most certainly not a Kooij one\n"
                         << "it has been corrected, but please, change that in your file.\n"
                         << "Thanks and a good day to you, user." << std::endl;
               kineticsModel = KineticsModel::ARRHENIUS;
            }
            if(verbose) 
              {
                std::cout << "   b: " << rate_constant->FirstChildElement("b")->GetText() << std::endl; 
              }
         }

//E has cal/mol default unit
          if(rate_constant->FirstChildElement("E") != NULL)
          {
            def_unit.set_unit("cal/mol");
            if(!rate_constant->FirstChildElement("E")->Attribute("units"))
            {
               antioch_unit_required("E",def_unit.get_symbol());
            }else
            {
               Units<NumericType> read_unit;
               read_unit.set_unit(rate_constant->FirstChildElement("E")->Attribute("units"));
               if(!read_unit.is_homogeneous(def_unit) &&
                  !read_unit.is_homogeneous("K")) //K directly given
               {
                   std::string errorstring("Error in reaction " + my_rxn->equation());
                   errorstring += "\n E units should be homogeneous to " + def_unit.get_symbol() + 
                                  " or K, and you provided " + read_unit.get_symbol();
                   antioch_unit_error(errorstring);
               }
               def_unit = read_unit;
            }
            data.push_back(std::atof(rate_constant->FirstChildElement("E")->GetText()));
            if(verbose) 
              {
                std::cout << "   E: " << rate_constant->FirstChildElement("E")->GetText()
                          << " " << def_unit.get_symbol() << std::endl; 
              }
          }


          if(rate_constant->FirstChildElement("D") != NULL)
          {
            def_unit.set_unit("K-1");
            if(!rate_constant->FirstChildElement("D")->Attribute("units"))
            {
               antioch_unit_required("D",def_unit.get_symbol());
            }else
            {
               Units<NumericType> read_unit;
               read_unit.set_unit(rate_constant->FirstChildElement("D")->Attribute("units"));
               if(!read_unit.is_homogeneous(def_unit))
               {
                 std::string errorstring("Error in reaction " + my_rxn->equation());
                 errorstring += "\n D units should be homogeneous to " + def_unit.get_symbol() + 
                                " and you provided " + read_unit.get_symbol();
                 antioch_unit_error(errorstring);
               }
               def_unit = read_unit;
            }
             data.push_back(std::atof(rate_constant->FirstChildElement("D")->GetText()) * def_unit.get_SI_factor());
            if(verbose) 
              {
                std::cout << "   D: " << rate_constant->FirstChildElement("D")->GetText()
                          << " " << def_unit.get_symbol() << std::endl; 
              }
          }

          //Tref
          if(kineticsModel == KineticsModel::HERCOURT_ESSEN ||
             kineticsModel == KineticsModel::BHE            ||
             kineticsModel == KineticsModel::KOOIJ          ||
             kineticsModel == KineticsModel::VANTHOFF) 
          {
            data.push_back(1.);
            def_unit.set_unit("K");
            if(!rate_constant->FirstChildElement("Tref"))
            {
                antioch_parameter_required("Tref","1 K");
            }else
            {
              if(!rate_constant->FirstChildElement("Tref")->Attribute("units"))
              {
                 antioch_unit_required("Tref",def_unit.get_symbol());
              }else
              {
                 Units<NumericType> read_unit;
                 read_unit.set_unit(rate_constant->FirstChildElement("Tref")->Attribute("units"));
                 if(!read_unit.is_homogeneous(def_unit))
                 {
                   std::string errorstring("Error in reaction " + my_rxn->equation());
                   errorstring += "\n Tref units should be homogeneous to " + def_unit.get_symbol() + 
                                  " and you provided " + read_unit.get_symbol();
                   antioch_unit_error(errorstring);
                }
                def_unit = read_unit;
              }
              data.back() = std::atof(rate_constant->FirstChildElement("Tref")->GetText()) * def_unit.get_SI_factor();
            }
          }

          //scale E -> E/R
          if(kineticsModel == KineticsModel::ARRHENIUS ||
             kineticsModel == KineticsModel::KOOIJ     ||
             kineticsModel == KineticsModel::VANTHOFF)
          {
            if(rate_constant->FirstChildElement("E")->Attribute("units"))
            {
               def_unit.set_unit("cal/mol"); // find E unit
               Units<NumericType> read_unit;
               read_unit.set_unit(rate_constant->FirstChildElement("E")->Attribute("units"));
               if(def_unit.is_homogeneous(read_unit)) //energy given
               {
                  def_unit.set_unit(read_unit.get_symbol() + "/K"); //unit of R is (E unit)/K
                  data.push_back(Antioch::Constants::R_universal<NumericType>() * Antioch::Constants::R_universal_unit<NumericType>().factor_to_some_unit(def_unit));
               }else if(read_unit.is_homogeneous("K")) //K directly given
               {
                  data.push_back(1.L);
               }else
               {
                   std::string errorstring("Error in reaction " + my_rxn->equation());
                   errorstring += "\n E units should be homogeneous to " + def_unit.get_symbol() + 
                                  " or K, and you provided " + read_unit.get_symbol();
                   antioch_unit_error(errorstring);
               }
            }else 
            {
               def_unit.set_unit("cal/mol/K"); // default R unit
               data.push_back(Antioch::Constants::R_universal<NumericType>() * Antioch::Constants::R_universal_unit<NumericType>().factor_to_some_unit(def_unit));
            }
          }

          //photochemistry
          // lambda is either a length (def nm) or cm-1
          // cross-section has several possibilities if given
          //   * cm2 per bin:
          //            - length (typically nm) or cm-1
          //   * cm2 no bin given: 
          //            - if given, lambda unit
          //            - if not, nm

          // starting with lambda (for bin unit in cross-section)
          // lambda is not in SI (m is really to violent), it will be nm
          if(rate_constant->FirstChildElement("lambda"))
          {
             data.clear();
             antioch_assert_equal_to(kineticsModel,KineticsModel::PHOTOCHEM);
             std::vector<std::string> lambda;
             def_unit.set_unit("nm");

       //reading part
            SplitString( std::string(rate_constant->FirstChildElement("lambda")->GetText()),
                          " ",
                          lambda,
                          /* include_empties = */ false );

         //unit checking part
             if(!rate_constant->FirstChildElement("lambda")->Attribute("units"))
             {
                antioch_unit_required("lambda",def_unit.get_symbol());
                for(unsigned int il = 0; il < lambda.size(); il++)
                {
                   data.push_back(std::atof(lambda[il].c_str()) * def_unit.factor_to_some_unit("nm"));
                }
             }else
             {
               Units<NumericType> read_unit;
               read_unit.set_unit(rate_constant->FirstChildElement("lambda")->Attribute("units"));
               if(read_unit.is_homogeneous(def_unit))
               {
                   def_unit = read_unit;
                   for(unsigned int il = 0; il < lambda.size(); il++)
                   {
                     data.push_back(std::atof(lambda[il].c_str()) * def_unit.factor_to_some_unit("nm"));
                   }
               }else if(read_unit.is_homogeneous("cm-1"))
               {
                 def_unit = read_unit;
                 for(unsigned int il = 0; il < lambda.size(); il++)
                 {
                    data.push_back(1.L/(std::atof(lambda[il].c_str()) * def_unit.factor_to_some_unit("nm-1")));
                  }
               }else
               {
                   std::string errorstring("Error in reaction " + my_rxn->equation());
                   errorstring += "\n Wavelength units should be homogeneous to " + def_unit.get_symbol() + " or cm-1"
                                  ", and you provided " + read_unit.get_symbol();
                   antioch_unit_error(errorstring);
               }
             }
             Antioch::Units<NumericType> bin_unit = def_unit;
             if(!rate_constant->FirstChildElement("cross_section"))
             {
                std::cerr << "Where is the cross-section?  In what universe have you photochemistry with a wavelength grid and no cross-section on it?" << std::endl;
                antioch_error();
             }

             /* here we will use two def unit:
              * cs_unit, cm2 by default
              * bin_unit, nm by default.
              *
              * strict rigorous unit is
              *         - (cs_unit - bin_unit): cm2/nm
              * correct unit is
              *         - cs_unit: cm2
              *
              * so we need to test against those two possibilities.
              * Now the funny part is that we test homogeneity, not
              * equality, for generality purposes, so in case of strict
              * rigorous unit, we need to decompose the read_unit into
              * cross_section and bin units, so we can make the appropriate change.
              * 
              * !TODO make the decomposition instead of strict equality
              */

             Antioch::Units<NumericType> cs_unit("cm2");
             if(!rate_constant->FirstChildElement("cross_section")->Attribute("units"))
             {
                antioch_unit_required("cross_section",(cs_unit - bin_unit).get_symbol());
             }else
             {
               Units<NumericType> read_unit;
               read_unit.set_unit(rate_constant->FirstChildElement("cross_section")->Attribute("units"));
               if(read_unit.is_homogeneous(cs_unit - bin_unit)) // here test the rigorous unit: cm2/nm 
               {
        // work to do here for decomposition   !!!!  supposes strict equality => cm2 per bin only
               }else if(read_unit.is_homogeneous(cs_unit)) // here test the almost rigorous unit cm2
               {
                  cs_unit = read_unit;
               }else //nothing can save you now...
               {
                   std::string errorstring("Error in reaction " + my_rxn->equation());
                   errorstring += "\n Cross-section units should be homogeneous to " + def_unit.get_symbol() + 
                                  ", and you provided " + read_unit.get_symbol();
                   antioch_unit_error(errorstring);
               }
             }
             std::vector<std::string> sigma;
             SplitString( std::string(rate_constant->FirstChildElement("cross_section")->GetText()),
                          " ",
                          sigma,
                          /* include_empties = */ false );
             if(sigma.size() != lambda.size())
             {
                std::cerr << "Your cross-section vector and your lambda vector don't have the same size!\n"
                          << "What am I supposed to do with that?"
                          << std::endl;
                antioch_error();
             }

             if(bin_unit.is_homogeneous("nm")) //nm
             {
              for(unsigned int ics = 0; ics < sigma.size(); ics++)
              {
                data.push_back(std::atof(sigma[ics].c_str()) * cs_unit.get_SI_factor() / (bin_unit.factor_to_some_unit("nm")));
              }
             }else if(bin_unit.is_homogeneous("cm-1"))//cm-1
             {   
               for(unsigned int ics = 0; ics < sigma.size(); ics++)
               {
                 data.push_back(std::atof(sigma[ics].c_str()) * cs_unit.get_SI_factor() / bin_unit.factor_to_some_unit("nm-1"));
               }
             }else //WHAT ?!!??
             {
                antioch_error();
             }
          } //end photochemistry

          if(data.empty()) //replace the old "if no A parameters" as A is not required anymore
          {
                std::cerr << "Somehow, I have a bad feeling about a chemical reaction without any data parameters...\n"
                          << "This is too sad, I give up...\n"
                          << "Please, check the reaction " << my_rxn->equation() << " before coming back to me." << std::endl;
                antioch_error(); //HEY!!!
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
        if(my_rxn->type() == ReactionType::LINDEMANN_FALLOFF ||
           my_rxn->type() == ReactionType::TROE_FALLOFF)
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
                    if( chem_mixture.species_name_map().count( pair.first ) )
                      {
                        my_rxn->set_efficiency( pair.first,
                                               chem_mixture.species_name_map().find( pair.first )->second,
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
           antioch_assert_equal_to (ReactionType::TROE_FALLOFF, my_rxn->type());
           FalloffReaction<NumericType,TroeFalloff<NumericType> > *my_fall_rxn =
                static_cast<FalloffReaction<NumericType,TroeFalloff<NumericType> > *> (my_rxn);

           Units<NumericType> def_unit;

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
           }else
           {
              def_unit.set_unit("K");
              if(!troe->FirstChildElement("T3")->Attribute("units"))
              {
                 antioch_unit_required("T3",def_unit.get_symbol());
              }else
              {
                 Units<NumericType> read_unit;
                 read_unit.set_unit(rate_constant->FirstChildElement("T3")->Attribute("units"));
                 if(!read_unit.is_homogeneous(def_unit))
                 {
                   std::string errorstring("Error in reaction " + my_rxn->equation());
                   errorstring += "\n Tref units should be homogeneous to " + def_unit.get_symbol() + 
                                  " and you provided " + read_unit.get_symbol();
                 }
                 def_unit = read_unit;
              }
           }
           my_fall_rxn->F().set_T3(std::atof(troe->FirstChildElement("T3")->GetText()) * def_unit.get_SI_factor());

           if(!troe->FirstChildElement("T1"))
           {
                std::cerr << "T* parameter of Troe falloff missing!" << std::endl;
                antioch_error();
           }else
           {
              def_unit.set_unit("K");
              if(!troe->FirstChildElement("T1")->Attribute("units"))
              {
                 antioch_unit_required("T1",def_unit.get_symbol());
              }else
              {
                 Units<NumericType> read_unit;
                 read_unit.set_unit(rate_constant->FirstChildElement("T1")->Attribute("units"));
                 if(!read_unit.is_homogeneous(def_unit))
                 {
                   std::string errorstring("Error in reaction " + my_rxn->equation());
                   errorstring += "\n Tref units should be homogeneous to " + def_unit.get_symbol() + 
                                  " and you provided " + read_unit.get_symbol();
                 }
                 def_unit = read_unit;
              }
           }
           my_fall_rxn->F().set_T1(std::atof(troe->FirstChildElement("T1")->GetText()) * def_unit.get_SI_factor());

           if(troe->FirstChildElement("T2"))//T2 is optional
           {
             def_unit.set_unit("K");
             if(!troe->FirstChildElement("T2")->Attribute("units"))
             {
                antioch_unit_required("T2",def_unit.get_symbol());
             }else
             {
                Units<NumericType> read_unit;
                read_unit.set_unit(rate_constant->FirstChildElement("T2")->Attribute("units"));
                if(!read_unit.is_homogeneous(def_unit))
                {
                  std::string errorstring("Error in reaction " + my_rxn->equation());
                  errorstring += "\n Tref units should be homogeneous to " + def_unit.get_symbol() + 
                                 " and you provided " + read_unit.get_symbol();
                }
               def_unit = read_unit;
             }
             my_fall_rxn->F().set_T2(std::atof(troe->FirstChildElement("T2")->GetText()) * def_unit.get_SI_factor());
           }
        }

        reaction_set.add_reaction(my_rxn);     

        if(verbose) std::cout << "\n\n";

        // Go to the next reaction
        reaction = reaction->NextSiblingElement("reaction");
      }

    return;
  }
 
} // end namespace Antioch

#endif // ANTIOCH_REACTION_SET_DATA_XML_H
