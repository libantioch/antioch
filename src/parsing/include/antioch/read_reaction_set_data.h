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

#ifndef ANTIOCH_READ_REACTION_SET_DATA_H
#define ANTIOCH_READ_REACTION_SET_DATA_H

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/string_utils.h"
#include "antioch/reaction_set.h"
#include "antioch/kinetics_parsing.h"
#include "antioch/reaction_parsing.h"
#include "antioch/physical_constants.h"
#include "antioch/units.h"
#include "antioch/parsing_enum.h"

// C++
#include <string>
#include <vector>
#include <map>

namespace Antioch
{

  template <typename NumericType>
  class ASCIIParser;

  template <typename NumericType>
  class ChemKinParser;

  template <typename NumericType>
  class ChemKinParser;

 /*!\file read_reaction_set_data.h
  *
  * We parse the file here, with an exhaustive
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

  template<class NumericType>
  void read_reaction_set_data_chemkin( const std::string& filename,
                                       const bool verbose,
                                       ReactionSet<NumericType>& reaction_set );


  template<typename NumericType>
  void read_reaction_set_data(const std::string &filename,
                              const bool verbose,
                              ReactionSet<NumericType>& reaction_set,
                              ParsingType type = ASCII );

  template <typename NumericType>
  void verify_unit_of_parameter(Units<NumericType> & default_unit, const std::string & provided_unit,
                                 const std::vector<std::string> & accepted_unit,
                                 const std::string & equation,     const std::string & parameter_name);


  /* ------------------------- Inline Functions -------------------------*/
  template<class NumericType>
  inline
  void read_reaction_set_data_xml( const std::string& filename,
                                   const bool verbose,
                                   ReactionSet<NumericType>& reaction_set )
  {
     read_reaction_set_data<NumericType >(filename,verbose,reaction_set,XML);
  }

  template<class NumericType>
  inline
  void read_reaction_set_data_chemkin( const std::string& filename,
                                       const bool verbose,
                                       ReactionSet<NumericType>& reaction_set )
  {
     read_reaction_set_data<NumericType>(filename,verbose,reaction_set,CHEMKIN);
  }

  template <typename NumericType>
  inline
  void verify_unit_of_parameter(Units<NumericType> & default_unit, const std::string & provided_unit,
                                const std::vector<std::string> & accepted_unit,
                                const std::string & equation, const std::string & parameter_name)
  {
     if(provided_unit.empty() && default_unit.is_united()) //no unit provided and there is one
     {
         antioch_unit_required(parameter_name,default_unit.get_symbol());
     }else // test
     {
        bool found(false);
        for(unsigned int i = 0; i < accepted_unit.size(); i++)
        {
           if(default_unit.is_homogeneous(provided_unit))
           {
              found = true;
              break;
           }
        }
        if(!found)
        {
           std::string errorstring("Error in reaction " + equation);
           errorstring += "\n" + parameter_name + " parameter's unit should be homogeneous to \"" + accepted_unit[0] + "\"";
           for(unsigned int i = 1; i < accepted_unit.size(); i++)errorstring += ", or \"" + accepted_unit[i] + "\"";
           errorstring += " and you provided \"" + provided_unit + "\"";
           antioch_unit_error(errorstring);
        }
        default_unit.set_unit(provided_unit);
     }
  }

  template<typename NumericType>
  inline
  void read_reaction_set_data( const std::string& filename,
                               const bool verbose,
                               ReactionSet<NumericType>& reaction_set,
                               ParsingType type )
  {
    ParserBase<NumericType> * parser(NULL);
    switch(type)
    {
      case ASCII:
         parser = new ASCIIParser<NumericType>(filename,verbose);
         break;
      case CHEMKIN:
         parser = new ChemKinParser<NumericType>(filename,verbose);
         break;
      case XML:
         parser = new XMLParser<NumericType>(filename,verbose);
         break;
      default:
         antioch_parsing_error("unknown type");
    }

    //error or no reaction data
    if(!parser->initialize())return;

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

    std::map<KineticsModel::KineticsModel,unsigned int> kinetics_model_map;
    kinetics_model_map[KineticsModel::CONSTANT]       = 0;
    kinetics_model_map[KineticsModel::HERCOURT_ESSEN] = 1;
    kinetics_model_map[KineticsModel::BERTHELOT]      = 2;
    kinetics_model_map[KineticsModel::ARRHENIUS]      = 3;
    kinetics_model_map[KineticsModel::BHE]            = 4;
    kinetics_model_map[KineticsModel::KOOIJ]          = 5;
    kinetics_model_map[KineticsModel::VANTHOFF]       = 7;
    kinetics_model_map[KineticsModel::PHOTOCHEM]      = 8;

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

    while (parser->reaction())
      {
        if (verbose) std::cout << "Reaction #" << parser->reaction_id() << ":\n"
                               << " eqn: " << parser->reaction_equation()
                               << std::endl;

        ReactionType::ReactionType typeReaction(ReactionType::ELEMENTARY);
        KineticsModel::KineticsModel kineticsModel(KineticsModel::HERCOURT_ESSEN); // = 0

        if (!parser->reaction_chemical_process().empty())
          {
            if (verbose) std::cout << " type: " << parser->reaction_chemical_process() << std::endl;
            if(!proc_keyword.count(parser->reaction_chemical_process()))
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
            typeReaction = proc_keyword[parser->reaction_chemical_process()];
          }
            
        bool reversible(parser->reaction_reversible());
        if (verbose) std::cout << "reversible: " << reversible << std::endl;
       
        kineticsModel = kin_keyword[parser->reaction_kinetics_model(models)];
        const std::string reading_kinetics_model = parser->reaction_kinetics_model(models);

        // usually Kooij is called Arrhenius, check here
        if(kineticsModel == KineticsModel::ARRHENIUS)
        {
          if(parser->verify_Kooij_in_place_of_Arrhenius())
          {
               kineticsModel = KineticsModel::KOOIJ;
               std::cerr << "In reaction " << parser->reaction_id() << "\n"
                         << "An equation of the form \"A * (T/Tref)^beta * exp(-Ea/(R*T))\" is a Kooij equation,\n"
                         << "I guess a modified Arrhenius could be a name too.  Whatever, the correct label is\n"
                         << "\"Kooij\", or, << à la limite >> \"ModifiedArrhenius\".  Please use those terms instead,\n"
                         << "thanks and a good day to you, user." << std::endl;
          }
        }

        // construct a Reaction object    
        Reaction<NumericType>* my_rxn = build_reaction<NumericType>(n_species, parser->reaction_equation(),
                                                                               reversible,typeReaction,kineticsModel);

        // We will add the reaction, unless we do not have a 
        // reactant or product
        bool relevant_reaction = true;
        unsigned int order_reaction(0);
        std::vector<std::pair<std::string,int> > molecules_pairs;

        if(parser->reactants_pairs(molecules_pairs))
          {
            if (verbose)
            {
              std::cout << "\n   reactants: "; 
              for(unsigned int ir = 0; ir < molecules_pairs.size(); ir++)
                  std::cout << molecules_pairs[ir].first << ":" << molecules_pairs[ir].second << " ";
            }
                
            for( unsigned int p=0; p < molecules_pairs.size(); p++ )
              {
                if(molecules_pairs[p].first == "e-") molecules_pairs[p].first = "e";

                if(verbose) std::cout  << "\n    " << molecules_pairs[p].first << " " << molecules_pairs[p].second;

                if( !chem_mixture.active_species_name_map().count( molecules_pairs[p].first ) )
                  {
                    relevant_reaction = false;
                    if (verbose) std::cout << "\n     -> skipping this reaction (no reactant " << molecules_pairs[p].first << ")";
                  }
                else
                  {
                    my_rxn->add_reactant( molecules_pairs[p].first,
                                         chem_mixture.active_species_name_map().find( molecules_pairs[p].first )->second,
                                         molecules_pairs[p].second );
                    order_reaction += molecules_pairs[p].second;
                  }
              }
          }

        molecules_pairs.clear();
        if(parser->products_pairs(molecules_pairs))
          {
            if(verbose) std::cout << "\n   products: ";
              for(unsigned int ir = 0; ir < molecules_pairs.size(); ir++)
                  std::cout << molecules_pairs[ir].first << ":" << molecules_pairs[ir].second << " ";
                
            for (unsigned int p=0; p < molecules_pairs.size(); p++)
              {
                if(molecules_pairs[p].first == "e-") molecules_pairs[p].first = "e";

                if(verbose) std::cout  << "\n    " << molecules_pairs[p].first << " " << molecules_pairs[p].second;

                if( !chem_mixture.active_species_name_map().count( molecules_pairs[p].first ) )
                  {
                    relevant_reaction = false;
                    if (verbose) std::cout << "\n     -> skipping this reaction (no product " << molecules_pairs[p].first << ")";
                  }
                else
                  {
                    my_rxn->add_product( molecules_pairs[p].first,
                                        chem_mixture.active_species_name_map().find( molecules_pairs[p].first )->second,
                                        molecules_pairs[p].second );
                  }
              }
             if(verbose) std::cout << std::endl;
          }

        if(!relevant_reaction)
        {
          if(verbose) std::cout << "skipped reaction\n\n";
          delete my_rxn;
          continue;
        }

        while(parser->rate_constant(reading_kinetics_model)) //for duplicate and falloff models, several kinetics rate to load, no mixing allowed
        {

         /* Any data is formatted by the parser method.
          * For any data required, parser sends back:
          *    - true/false if data is found
          *    - value of data
          *    - unit of data if found, empty string else
          *    - default unit of data
          *
          * The parser defines the defaults, note the special case
          * of the pre-exponential parameters:
          *  its unit is [quantity-1]^(order - 1)/s, thus the parser
          *  defines only the [quantity-1] unit (SI unit is m^3/mol)
          *  as default.
          */

          std::vector<NumericType> data; // for rate constant

          Units<NumericType> def_unit;
          int pow_unit(order_reaction - 1);
          //threebody always
          if(my_rxn->type() == ReactionType::THREE_BODY)pow_unit++;
          //falloff for k0
          if(my_rxn->type() == ReactionType::LINDEMANN_FALLOFF ||
             my_rxn->type() == ReactionType::TROE_FALLOFF)
          {
             //k0 is either determined by an explicit name, or is the first of unnamed rate constants
             if(parser->is_k0(my_rxn->n_rate_constants(),reading_kinetics_model))pow_unit++;
          }

          NumericType par_value(-1.);
          std::vector<NumericType> par_values;
          std::string par_unit;
          std::string default_unit;
          std::vector<std::string> accepted_unit;

        // verbose as we read along
          if(verbose)std::cout << " rate: " << models[kinetics_model_map[kineticsModel]] << " model\n";

// pre-exponential
          if(parser->rate_constant_preexponential_parameter(par_value, par_unit, default_unit))
          {
            // using Units object to build accepted_unit
            accepted_unit.clear();
            def_unit.set_unit("m3/mol");
   //to the m-1 power
            if(pow_unit != 0)
            {
               def_unit *= pow_unit;
            }else
            {
              def_unit.clear();
            }
            def_unit.substract("s");    // per second
            accepted_unit.push_back(def_unit.get_symbol());

            def_unit.set_unit(default_unit);
   //to the m-1 power
            if(pow_unit != 0)
            {
               def_unit *= pow_unit;
            }else
            {
              def_unit.clear();
            }
            def_unit.substract("s");    // per second
            verify_unit_of_parameter(def_unit, par_unit, accepted_unit, my_rxn->equation(), "A");
           if(verbose) 
             {
             std::cout  << "   A: " << par_value
                        << " "      << def_unit.get_symbol() << std::endl; 
             }
           data.push_back(par_value * def_unit.get_SI_factor());
         }

// beta
         if(parser->rate_constant_power_parameter(par_value,par_unit,default_unit))
         {
            accepted_unit.clear();
            accepted_unit.push_back("");
            def_unit.set_unit(default_unit);
            verify_unit_of_parameter(def_unit, par_unit, accepted_unit, my_rxn->equation(), "beta");
            if(par_value == 0.)//if ARRHENIUS parameterized as KOOIJ, bad test, need to rethink it
            {
               std::cerr << "In reaction " << parser->reaction_id() << "\n"
                         << "An equation of the form \"A * exp(-Ea/(R*T))\" is an Arrhenius equation,\n"
                         << "and most certainly not a Kooij one\n"
                         << "it has been corrected, but please, change that in your file.\n"
                         << "Thanks and a good day to you, user." << std::endl;
               kineticsModel = KineticsModel::ARRHENIUS;
            }else
            {
              data.push_back(par_value * def_unit.get_SI_factor());
            }
            if(verbose) 
              {
                std::cout << "   b: " << par_value << std::endl; 
              }
         }

// activation energy
          if(parser->rate_constant_activation_energy_parameter(par_value,par_unit,default_unit))
          {
            accepted_unit.clear();
            accepted_unit.push_back("J/mol");
            accepted_unit.push_back("K");
            def_unit.set_unit(default_unit);
            verify_unit_of_parameter(def_unit, par_unit, accepted_unit, my_rxn->equation(), "Ea");
            data.push_back(par_value * def_unit.get_SI_factor());
            if(verbose) 
              {
                std::cout << "   E: " << par_value
                          << " "      << def_unit.get_symbol() << std::endl; 
              }
          }


// Berthelot coefficient (D)
          if(parser->rate_constant_Berthelot_coefficient_parameter(par_value,par_unit,default_unit))
          {
            accepted_unit.clear();
            accepted_unit.push_back("K");
            def_unit.set_unit(default_unit);
            verify_unit_of_parameter(def_unit, par_unit, accepted_unit, my_rxn->equation(), "D");
            data.push_back(par_value * def_unit.get_SI_factor());
            if(verbose) 
              {
                std::cout << "   D: " << par_value
                          << " "      << def_unit.get_symbol() << std::endl; 
              }
          }

// Tref, not for everyone
          if(kineticsModel == KineticsModel::HERCOURT_ESSEN ||
             kineticsModel == KineticsModel::BHE            ||
             kineticsModel == KineticsModel::KOOIJ          ||
             kineticsModel == KineticsModel::VANTHOFF) 
          {
            par_value = 1.;
            if(parser->rate_constant_Tref_parameter(par_value,par_unit,default_unit))
            {
              accepted_unit.clear();
              accepted_unit.push_back("K");
              def_unit.set_unit(default_unit);
              verify_unit_of_parameter(def_unit, par_unit, accepted_unit, my_rxn->equation(), "Tref");
            }else
            {
                antioch_parameter_required("Tref","1 K");
            }
            data.push_back(par_value);
          }

// scale E -> E/R
          if(kineticsModel == KineticsModel::ARRHENIUS ||
             kineticsModel == KineticsModel::KOOIJ     ||
             kineticsModel == KineticsModel::VANTHOFF)
          {
            parser->rate_constant_activation_energy_parameter(par_value,par_unit,default_unit);
            (par_unit.empty())?def_unit.set_unit(default_unit):def_unit.set_unit(par_unit);
        // now finding R unit: [Ea] / [K]
            def_unit.substract("K");

            par_value = (def_unit.is_united())?
                                Antioch::Constants::R_universal<NumericType>() // Ea already tranformed in SI
                                              :1.L;  // no unit, so Ea already in K
            data.push_back(par_value);
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
          if(parser->rate_constant_lambda_parameter(par_values,par_unit,default_unit))
          {
             antioch_assert_equal_to(kineticsModel,KineticsModel::PHOTOCHEM);
             accepted_unit.clear();
             accepted_unit.push_back("nm");
             accepted_unit.push_back("cm-1");
             def_unit.set_unit(default_unit);
             verify_unit_of_parameter(def_unit, par_unit, accepted_unit, my_rxn->equation(), "lambda");
             
             data.clear();
         
             if(def_unit.is_homogeneous("nm"))// it's a length
             {
                for(unsigned int il = 0; il < par_values.size(); il++)
                {
                   data.push_back(par_values[il] * def_unit.factor_to_some_unit("nm"));
                }
             }else // it's the inverse of a length
             {
                for(unsigned int il = 0; il < par_values.size(); il++)
                {
                   data.push_back(1.L/(par_values[il] * def_unit.factor_to_some_unit("nm-1")));
                }
             }

             //now the cross-section

             NumericType bin_coefficient = (def_unit.is_homogeneous("nm"))?def_unit.factor_to_some_unit("nm"):
                                                                           def_unit.factor_to_some_unit("nm-1");
             if(!parser->rate_constant_cross_section_parameter(par_values,par_unit,default_unit))
             {
                std::cerr << "Where is the cross-section?  In what universe have you photochemistry with a wavelength grid and no cross-section on it?" << std::endl;
                antioch_error();
             }
             //test length
             if(par_values.size() != data.size())
             {
                std::cerr << "Your cross-section vector and your lambda vector don't have the same size!\n"
                          << "What am I supposed to do with that?"
                          << std::endl;
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

             accepted_unit.clear();
             accepted_unit.push_back("cm2");       // only cross-section
             accepted_unit.push_back("cm2/nm");    // per bin, bin is length-like
             accepted_unit.push_back("cm2/nm-1");  // per bin, bin is inverse length-like
             verify_unit_of_parameter(def_unit, par_unit, accepted_unit, my_rxn->equation(), "cross-section");
             

             if(def_unit.is_homogeneous("cm2")) // gotta use the bin unit, real unit is [cs]/[provided bin unit]
             {
               for(unsigned int ics = 0; ics < par_values.size(); ics++)
               {
                  data.push_back(par_values[ics] * def_unit.get_SI_factor() / bin_coefficient); //cs in SI, bin in nm or nm-1
               }
             }else                              // bin unit is provided
             {
               std::string target_unit = (def_unit.is_homogeneous("cm2/nm"))?"m2/nm":"m2/nm-1";
               for(unsigned int ics = 0; ics < par_values.size(); ics++)
               {
                  data.push_back(par_values[ics] * def_unit.factor_to_some_unit(target_unit));
               }
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
           if(parser->where_is_k0(reading_kinetics_model) == 1) // second given is k0
           { 
               my_rxn->swap_forward_rates(0,1);
           }
        }


        std::vector<std::pair<std::string,NumericType> > efficiencies;
        //efficiencies are only for three body reactions
        if(parser->efficiencies(efficiencies))
          {
             antioch_assert_equal_to (ReactionType::THREE_BODY, my_rxn->type());

             for(unsigned int p = 0; p < efficiencies.size(); p++)
               {
                    if(verbose)std::cout  << "\n" << efficiencies[p].first << " " << efficiencies[p].second;

                    if(efficiencies[p].first == "e-") efficiencies[p].first = "e";

                    // it is possible that the efficiency is specified for a species we are not
                    // modeling - so only add the efficiency if it is included in our list
                    if( chem_mixture.active_species_name_map().count( efficiencies[p].first ) )
                      {
                        my_rxn->set_efficiency( efficiencies[p].first,
                                               chem_mixture.active_species_name_map().find( efficiencies[p].first )->second,
                                               efficiencies[p].second );
                      }
               }
               if(verbose)std::cout << std::endl;
          }

        //F parameters only for Troe falloff
        if(parser->Troe())
        {
           antioch_assert_equal_to (ReactionType::TROE_FALLOFF, my_rxn->type());
           FalloffReaction<NumericType,TroeFalloff<NumericType> > *my_fall_rxn =
                static_cast<FalloffReaction<NumericType,TroeFalloff<NumericType> > *> (my_rxn);

           Units<NumericType> def_unit;
           NumericType par_value(-1.);
           std::string par_unit;
           std::string default_unit;
           std::vector<std::string> accepted_unit;

        // alpha
           if(!parser->Troe_alpha_parameter(par_value,par_unit,default_unit))
           {
                std::cerr << "alpha parameter of Troe falloff missing!" << std::endl;
                antioch_error();
           }
           accepted_unit.clear();
           accepted_unit.push_back("");
           def_unit.set_unit(default_unit);
           verify_unit_of_parameter(def_unit, par_unit, accepted_unit, my_fall_rxn->equation(), "alpha");
           my_fall_rxn->F().set_alpha(par_value * def_unit.get_SI_factor());

        // T***
           if(!parser->Troe_T3_parameter(par_value,par_unit,default_unit))
           {
                std::cerr << "T*** parameter of Troe falloff missing!" << std::endl;
                antioch_error();
           }
           accepted_unit.clear();
           accepted_unit.push_back("K");
           def_unit.set_unit(default_unit);
           verify_unit_of_parameter(def_unit, par_unit, accepted_unit, my_fall_rxn->equation(), "T***");
           my_fall_rxn->F().set_T3(par_value * def_unit.get_SI_factor());

        // T*
           if(!parser->Troe_T1_parameter(par_value,par_unit,default_unit))
           {
                std::cerr << "T* parameter of Troe falloff missing!" << std::endl;
                antioch_error();
           }
        // accepted unit is the same
           def_unit.set_unit(default_unit);
           verify_unit_of_parameter(def_unit, par_unit, accepted_unit, my_fall_rxn->equation(), "T*");
           my_fall_rxn->F().set_T1(par_value * def_unit.get_SI_factor());

        // T** is optional
           if(parser->Troe_T2_parameter(par_value,par_unit,default_unit))
           {
             def_unit.set_unit(default_unit);
             verify_unit_of_parameter(def_unit, par_unit, accepted_unit, my_fall_rxn->equation(), "T**");
             my_fall_rxn->F().set_T2(par_value * def_unit.get_SI_factor());
           }
        }

        reaction_set.add_reaction(my_rxn);     

        if(verbose) std::cout << "\n\n";
      }

    if(parser)delete parser;
    return;
  }
 
} // end namespace Antioch

#endif // ANTIOCH_READ_REACTION_SET_DATA_H
