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
#ifndef ANTIOCH_CHEMKIN_PARSER_H
#define ANTIOCH_CHEMKIN_PARSER_H

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/string_utils.h"
#include "antioch/parsing_enum.h"

//ChemKin

//C++
#include <sstream>
#include <string>
#include <vector>
#include <map>

namespace Antioch{


  /*! ChemKin format file reader
   *
   * defaults unit:
   *  - A: cm-molecule or cm-mol (see read_reaction_set_data.h for addition of s-1)
   *  - beta: no units
   *  - Ea: cal/mol
   *
   * There are no other kinetics paramters for rate constant. Chemical processes
   * are:
   *  - Elementary
   *  - Duplicate
   *  - Three-Body
   *  - Lindemann falloff
   *  - Troe falloff
   *  - SRI falloff (not supported)
   */
  template <typename NumericType = double>
  class ChemKinParser{
        public:
          ChemKinParser(const std::string &filename);
          ~ChemKinParser();

//// first local pointers
         /*! Read header of file, go to interesting part*/
         bool initialize();

         /*! go to next reaction*/
         bool reaction();

         /*! go to next rate constant*/
         bool rate_constant(const std::string & /* kinetics_model */);

         /*! return true if there's a Troe block*/
         bool Troe();

         /*! return reaction id, 0 if not provided*/
         const std::string reaction_id() const;

         /*! return reaction equation */
         const std::string reaction_equation() const;

         /*! return reaction chemical process*/
         const std::string reaction_chemical_process() const;

         /*! return reaction kinetics model*/
         const std::string reaction_kinetics_model(const std::vector<std::string> & /*kinetics_models*/) const;

         /*! return reversible state*/
         bool reaction_reversible() const;

         /*! return pairs of reactants and stoichiometric coefficients*/
         bool reactants_pairs(std::vector<std::pair<std::string,int> > & reactants_pair) const;

         /*! return pairs of products and stoichiometric coefficients*/
         bool products_pairs(std::vector<std::pair<std::string,int> > & products_pair) const;

         /*! return true if "name" attribute is found with value "k0"*/
         bool is_k0(unsigned int nrc, const std::string & kin_model) const;

         /*! return index of k0 (0 or 1)*/
         unsigned int where_is_k0(const std::string & /*kin_model*/) const;

         /*! return true if pre exponentiel coefficient*/
         bool rate_constant_preexponential_parameter(    NumericType & A,    std::string & A_unit,    std::string & def_unit) const;

         /*! return true if beta coefficient*/
         bool rate_constant_power_parameter(             NumericType & b,    std::string & b_unit,    std::string & def_unit) const;

         /*! return true if activation energie*/
         bool rate_constant_activation_energy_parameter(NumericType & Ea,   std::string & Ea_unit,   std::string & def_unit) const;

         /*! return true if D coefficient*/
         bool rate_constant_Berthelot_coefficient_parameter(NumericType & D,    std::string & D_unit,    std::string & def_unit) const;

         /*! return true if Tref*/
         bool rate_constant_Tref_parameter(              NumericType & Tref, std::string & Tref_unit, std::string & def_unit) const;

         /*! return true if lambda*/
         bool rate_constant_lambda_parameter(       std::vector<NumericType> & lambda, std::string & lambda_unit, std::string & def_unit) const;

         /*! return true if sigma*/
         bool rate_constant_cross_section_parameter(std::vector<NumericType> & sigma,  std::string & sigma_unit,  std::string & def_unit) const;

         /*! return true if a Kooij is called Arrhenuis*/
         bool verify_Kooij_in_place_of_Arrhenius() const;

         /*! return true if efficiencies are found*/
         bool efficiencies(std::vector<std::pair<std::string,NumericType> > & par_values) const;

         /*! return true is alpha*/
         bool Troe_alpha_parameter(NumericType & alpha, std::string & alpha_unit, std::string & def_unit) const;

         /*! return true is alpha*/
         bool Troe_T1_parameter(   NumericType & T1,    std::string & T1_unit,    std::string & def_unit) const;

         /*! return true is alpha*/
         bool Troe_T2_parameter(   NumericType & T2,    std::string & T2_unit,    std::string & def_unit) const;

         /*! return true is alpha*/
         bool Troe_T3_parameter(   NumericType & T3,    std::string & T3_unit,    std::string & def_unit) const;

        private:

          /*! Convenient method */
          void parse_a_line(const std::string & line);

          /*! Convenient method */
          void parse_equation_coef(const std::string & line);

          /*! Convenient method */
          void parse_coefficients_line(const std::string &line);

          /*! Convenient method */
          std::pair<std::string,NumericType> parse_molecule(const std::string & molecule);

          /*! Convenient method */
          bool is_real_number(const char & c) const;

          /*! if stoichiometry is real, make them integer*/
          void rescale_stoichiometry();

          /*! check if the stoichiometric coef is an integer */
          bool after_coma_digits(NumericType number) const;

          /*! look for q when given r = p / q*/
          unsigned int factor_to_int(NumericType number) const;

          /*! Never use default constructor*/
          ChemKinParser();
          std::ifstream                    _doc;

          bool                             _reversible;
          unsigned int                     _nrates; // total number of rates
          unsigned int                     _crates; // current place

          // ChemKin allows real stoichiometric coefficients
          std::vector<std::pair<std::string,NumericType> > _reactants;
          std::vector<std::pair<std::string,NumericType> > _products;

          std::string                      _equation;
          std::string                      _chemical_process;
          std::string                      _kinetics_model;

          std::vector<NumericType>         _A;
          std::vector<NumericType>         _b;
          std::vector<NumericType>         _Ea;
          std::vector<NumericType>         _D;     // ChemKin don't use it, but whatever, let's generalize, just in case

          std::vector<std::pair<std::string,NumericType> > _efficiencies;

          NumericType                      _Tref;  // ChemKin don't use it, but whatever, let's generalize, just in case
          NumericType                      _Troe_alpha;
          NumericType                      _Troe_T1;
          NumericType                      _Troe_T2;
          NumericType                      _Troe_T3;

          std::map<ParsingKey,std::string> _map;
          std::map<ParsingKey,std::string> _default_unit;

          std::map<std::string,std::string> _unit_custom_ea;
          std::map<std::string,std::string> _unit_custom_A;

/*ChemKin spec*/

          class ChemKinSpec
          {
            public:
              ChemKinSpec():_duplicate("DUPLICATE"),_end_tag("END"),_comment("!"),_parser("/")
                           {
                               _delim[PLUS]           = "+";
                               _delim[REVERSIBLE]     = "=";
                               _delim[REVERSIBLE_ALT] = "<=>";
                               _delim[IRREVERSIBLE]   = "=>";

                               _symbol[TB]      = "+M";
                               _symbol[FAL]     = "(+M)";
                               _symbol[PHOTO]   = "HV";
                               _symbol[ELECTRO] = "E";
                           };

              ~ChemKinSpec(){}

              enum Delim{
                          PLUS,
                          IRREVERSIBLE,
                          REVERSIBLE,
                          REVERSIBLE_ALT
                        };

              enum Symbol{
                          TB,
                          FAL,
                          PHOTO,
                          ELECTRO 
                         };

              const std::map<Delim,std::string> & delim()   const {return _delim;}

              const std::map<Symbol,std::string> & symbol() const {return _symbol;}

              const std::string & duplicate()               const {return _duplicate;}

              const std::string & end_tag()                 const {return _end_tag;}

              const std::string & comment()                 const {return _comment;}

              const std::string & parser()                  const {return _parser;}

              bool is_comment(const char & c)               const {return (c == _comment[0]);}

              bool is_equation_delimiter(const std::string & test) const {return (test == _delim.at(REVERSIBLE)     || 
                                                                                  test == _delim.at(REVERSIBLE_ALT) || 
                                                                                  test == _delim.at(IRREVERSIBLE)
                                                                                  );}

            private:
              const std::string            _duplicate;
              const std::string            _end_tag;
              const std::string            _comment;
              const std::string            _parser;
              std::map<Delim,std::string>  _delim;
              std::map<Symbol,std::string> _symbol;

          };

          const ChemKinSpec _spec;

  };


  template <typename NumericType>
  inline
  ChemKinParser<NumericType>::ChemKinParser(const std::string &filename)
  {
     _doc.open(filename.c_str());
    if(_doc.fail())
      {
        std::cerr << "ERROR: unable to load ChemKin file " << filename << std::endl;
        antioch_error();
      }

      _map[ParsingKey::REACTION_DATA]    = "REAC"; //REACTIONS || REAC
      _map[ParsingKey::FALLOFF_LOW_NAME] = "LOW";
      _map[ParsingKey::TROE_FALLOFF]     = "TROE";

   // typically chemkin files list 
   //      pre-exponential parameters in (m3/kmol)^(m-1)/s
   //      activation energy in cal/mol, but we want it in K.
   //      power parameter without unit
   // if falloff, we need to know who's k0 and kinfty
   // if photochemistry, we have a cross-section on a lambda grid
   //      cross-section typically in cm2/nm (cross-section on a resolution bin, 
   //                                          if bin unit not given, it is lambda unit (supposed to anyway), and a warning message)
   //      lambda typically in nm, sometimes in ang, default considered here is nm
   //                         you can also have cm-1, conversion is done with
   //                         formulae nm = cm-1 * / * adapted factor
      _default_unit[ParsingKey::PREEXP]                = "cm3/mol";
      _default_unit[ParsingKey::POWER]                 = "";
      _default_unit[ParsingKey::ACTIVATION_ENERGY]     = "cal/mol";
      _default_unit[ParsingKey::BERTHELOT_COEFFICIENT] = "K-1";
      _default_unit[ParsingKey::TREF]                  = "K";
      _default_unit[ParsingKey::HV_LAMBDA]             = "nm";
      _default_unit[ParsingKey::HV_CROSS_SECTION]      = "cm2/nm";
      _default_unit[ParsingKey::EFFICIENCY]            = "";
      _default_unit[ParsingKey::TROE_F_ALPHA]          = "";
      _default_unit[ParsingKey::TROE_F_TS]             = "K";
      _default_unit[ParsingKey::TROE_F_TSS]            = "K";
      _default_unit[ParsingKey::TROE_F_TSSS]           = "K";

      _unit_custom_ea["CAL/MOL"]                       = "cal/mol";
      _unit_custom_ea["KCAL/MOL"]                      = "kcal/mol";
      _unit_custom_ea["JOULES/MOL"]                    = "J/mol";
      _unit_custom_ea["KELVINS"]                       = "K";
      _unit_custom_ea["EVOLTS"]                        = "eV";

      _unit_custom_A["MOLES"]                          = "cm3/mol";
      _unit_custom_A["MOLECULES"]                      = "cm3/molecule";

  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::initialize()
  {
      std::string line;
      std::getline(_doc,line);
      bool init(true);

      while(line.find(_map.at(ParsingKey::REACTION_DATA)) == std::string::npos)
      {
         if(!std::getline(_doc,line) || _doc.eof())
         {
            init = false;
            break;
         }
      }

      if(init)
      {
         std::vector<std::string> keywords;
         int nw = SplitString(line," ",keywords,false);
         if(nw < 1)antioch_parsing_error("ChemKin parser: now you're in trouble, can find the REACTION tag anymore.");
         for(unsigned int k = 1; k < keywords.size(); k++)// first word is REACTION
         {
            if(_unit_custom_ea.count(keywords[k]))
            {
               _default_unit[ParsingKey::PREEXP]  = _unit_custom_ea.at(keywords[k]);
            }else if(_unit_custom_A.count(keywords[k]))
            {
               _default_unit[ParsingKey::ACTIVATION_ENERGY] = _unit_custom_A.at(keywords[k]);
            }else
            {
               antioch_parsing_error("ChemKin parser: I don't have ChemKin supporting this word as a unit:\n" + keywords[k]);
            }
         }
      }

      return init;
  }

  template <typename NumericType>
  inline
  ChemKinParser<NumericType>::~ChemKinParser()
  {
     _doc.close();
     return;
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::reaction()
  {
      /*default*/
      _reversible = true;
      _nrates     = 0;
      _crates     = 0;

      _kinetics_model   = "Kooij"; //always Kooij, sometimes Arrhenius, should be corrected within the parser
      _chemical_process = "Elementary";
      _equation.clear();

      _reactants.clear();
      _products.clear();

      _A.clear();
      _b.clear();
      _Ea.clear();
      _D.clear();

      _Tref       =  1.;
      _Troe_alpha = -1.;
      _Troe_T1    = -1.;
      _Troe_T2    = -1.;
      _Troe_T3    = -1.;

      /* reaction */
      bool reac = true;
      std::string line;
      std::getline(_doc,line);
      while(!line.empty())
      {
        if(_spec.is_comment(line[0]))continue;
        if(line.find(_spec.end_tag()) != std::string::npos || _doc.eof()) // equivalent
        {
           reac = false;
           break;
        }
        if(line.find(_spec.comment()) != std::string::npos)line.erase(line.find(_spec.comment()),std::string::npos); //supress comment
        this->parse_a_line(line);
        std::getline(_doc,line);
      }

      return reac;
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::rate_constant(const std::string & /* kinetics_model */)
  {
     _crates++;
     return (_crates <= _nrates);
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::Troe()
  {
      return (_chemical_process == "TroeFalloff");
  }

  template <typename NumericType>
  inline
  const std::string ChemKinParser<NumericType>::reaction_id() const
  {
     return _equation;
  }

  template <typename NumericType>
  inline
  const std::string ChemKinParser<NumericType>::reaction_equation() const
  {
     return _equation;
  }

  template <typename NumericType>
  inline
  const std::string ChemKinParser<NumericType>::reaction_chemical_process() const
  {
      return _chemical_process;
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::reaction_reversible() const
  {
     return _reversible;
  }

  template <typename NumericType>
  inline
  const std::string ChemKinParser<NumericType>::reaction_kinetics_model(const std::vector<std::string> & /*kinetics_models*/) const
  {
      return _kinetics_model;
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::reactants_pairs(std::vector<std::pair<std::string,int> > & reactants_pair) const
  {
      reactants_pair.clear();
      reactants_pair.resize(_reactants.size());
      for(unsigned int i = 0; i < _reactants.size(); i++)
      {
        reactants_pair[i] = std::make_pair(_reactants[i].first,(int)_reactants[i].second);
      }
      return !(_reactants.empty());
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::products_pairs(std::vector<std::pair<std::string,int> > & products_pair) const
  {
      products_pair.clear();
      products_pair.resize(_products.size());
      for(unsigned int i = 0; i < _products.size(); i++)
      {
        products_pair[i] = std::make_pair(_products[i].first,(int)_products[i].second);
      }
      return !(_products.empty());
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::is_k0(unsigned int nrc, const std::string & kin_model) const
  {
// k0 is always the first, explicit in ChemKin
      return (nrc == 0);
  }

  template <typename NumericType>
  inline
  unsigned int ChemKinParser<NumericType>::where_is_k0(const std::string & /*kin_model*/) const
  {
     return 0;
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::rate_constant_preexponential_parameter(NumericType & A, std::string & A_unit, std::string & def_unit) const
  {
     if(_crates <= _A.size())
     {
       A = _A[_crates - 1];
// there is no default unit (or always default), anyway, units are
// always explicit
       A_unit = _default_unit.at(ParsingKey::PREEXP);
       def_unit = A_unit;
     }
     return (_crates <= _A.size());
  }


  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::rate_constant_power_parameter(NumericType & b, std::string & b_unit, std::string & def_unit) const
  {
    if(_crates <= _b.size())
    {
       b = _b[_crates-1];
// there is no default unit (or always default), anyway, units are
// always explicit
       b_unit = _default_unit.at(ParsingKey::POWER);
       def_unit = b_unit;
     }
     return (_crates <= _b.size());
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::rate_constant_activation_energy_parameter(NumericType & Ea, std::string & Ea_unit, std::string & def_unit) const
  {
     if(_crates <= _Ea.size());
     {
       Ea = _Ea[_crates-1];
// there is no default unit (or always default), anyway, units are
// always explicit
       Ea_unit = _default_unit.at(ParsingKey::ACTIVATION_ENERGY);
       def_unit = Ea_unit;
     }
     return (_crates <= _Ea.size());
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::rate_constant_Berthelot_coefficient_parameter(NumericType & D, std::string & D_unit, std::string & def_unit) const
  {
      return false; //not a chemkin model
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::rate_constant_lambda_parameter(std::vector<NumericType> & lambda, std::string & lambda_unit, std::string & def_unit) const
  {
      return false; //not a supported chemkin model yet
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::rate_constant_cross_section_parameter(std::vector<NumericType> & sigma, std::string & sigma_unit, std::string & def_unit) const
  {
      return false; //not a supported chemkin model yet
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::rate_constant_Tref_parameter(NumericType & Tref, std::string & Tref_unit, std::string & def_unit) const
  {
     Tref = 1.L;
     Tref_unit = _default_unit.at(ParsingKey::TREF);
     def_unit = Tref_unit;
     return (_crates <= _b.size());
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::verify_Kooij_in_place_of_Arrhenius() const
  {
      return false;
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::efficiencies(std::vector<std::pair<std::string,NumericType> > & par_values) const
  {
     par_values = _efficiencies;
     return !_efficiencies.empty();
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::Troe_alpha_parameter(NumericType & alpha, std::string & alpha_unit, std::string & def_unit) const
  {
    alpha = _Troe_alpha;
    alpha_unit = _default_unit.at(ParsingKey::TROE_F_ALPHA);
    def_unit = alpha_unit;
    
    return (_chemical_process == "TroeFalloff");
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::Troe_T1_parameter(NumericType & T1, std::string & T1_unit, std::string & def_unit) const
  {
    T1 = _Troe_T1;
    T1_unit = _default_unit.at(ParsingKey::TROE_F_TS);
    def_unit = T1_unit;

    return (_chemical_process == "TroeFalloff");
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::Troe_T2_parameter(NumericType & T2, std::string & T2_unit, std::string & def_unit) const
  {
    T2 = _Troe_T2;
    T2_unit = _default_unit.at(ParsingKey::TROE_F_TSS);
    def_unit = T2_unit;

    return (_Troe_T2 > 0.);
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::Troe_T3_parameter(NumericType & T3, std::string & T3_unit, std::string & def_unit) const
  {
    T3 = _Troe_T3;
    T3_unit = _default_unit.at(ParsingKey::TROE_F_TSSS);
    def_unit = T3_unit;

    return (_chemical_process == "TroeFalloff");
  }

  template <typename NumericType>
  inline
  void ChemKinParser<NumericType>::parse_a_line(const std::string & line)
  {
     if(line.find(_spec.delim().at(_spec.REVERSIBLE)) != std::string::npos) //equation a beta ea
     {
        this->parse_equation_coef(line);
     }
     else if(line.find(_spec.duplicate()) != std::string::npos)
     {
        _chemical_process = "Duplicate";
     }else
     {
        this->parse_coefficients_line(line);
     }
  }


  template <typename NumericType>
  inline
  void ChemKinParser<NumericType>::parse_equation_coef(const std::string & line)
  {

     std::vector<std::string> out;
     int nword = SplitString(line," ",out,false);
     if(nword < 4)antioch_parsing_error("ChemKin parser: unrecognized reaction input line:\n" + line);

// parameters are the three last components
     _Ea.push_back(std::atof(out[nword-1].c_str())); // last
     _b.push_back( std::atof(out[nword-2].c_str())); // last - 1 
     _A.push_back( std::atof(out[nword-3].c_str())); // last - 2

// equation is the rest, can be 1 word as nreac + nprod + 1
// so making it 1 whatever happens
     std::string equation;
     for(unsigned int i = 0; i < nword - 3; i++)
     {
        equation += out[i];
     }

// first we need to treat the (+M) case (falloff)
// as it is not compatible with SplitString using ChemKinSpec::PLUS (+)
// as delimiter
     if(equation.find(_spec.symbol().at(ChemKinSpec::FAL)) != std::string::npos)
     {
        // Lindemann by default
        _chemical_process = "LindemannFalloff";
        // supress all occurences
        while(equation.find(_spec.symbol().at(ChemKinSpec::FAL)) != std::string::npos)
        {
           equation.erase(equation.find(_spec.symbol().at(ChemKinSpec::FAL)),4);
        }
     }

     out.clear();
     nword = SplitString(line,_spec.delim().at(ChemKinSpec::PLUS),out,true); //empties are cations charge, formatted as reac_i equation_separator prod_i
     if(nword < 3)antioch_parsing_error("ChemKin parser: unrecognized reaction equation:\n" + equation);
/* cases are:
    - equation_separator => switch between reac and prod 
    - molecules:
        * reactant or prod, cation if followed by empty
        * M if three body
*/
     bool prod(false);
     for(unsigned int i = 0; i < nword; i++)
     {
        if(_spec.is_equation_delimiter(out[i]))
        {
           _reversible = !(out[i] == _spec.delim().at(ChemKinSpec::IRREVERSIBLE));
           prod = true;
         // is it a third-body reaction?
        }else if(out[i] == _spec.symbol().at(ChemKinSpec::TB))
        {
           if(_chemical_process.find("Falloff") != std::string::npos)
                      antioch_parsing_error("ChemKin parser: it seems you want both a falloff and a three-body reaction:\n" + equation);

           _chemical_process = "ThreeBoby";

        }else // here's a regular molecule
        {
         // no assumption on the charge, you can put as many '+' as you want
         // adding empties, then skipping them
           unsigned int j(1);
           while(out[i+j].empty())
           {
             out[i] += "+";
             j++;
           }
           j--; // back to non empty
           (prod)?_products.push_back(this->parse_molecule(out[i]))
                  :
                  _reactants.push_back(this->parse_molecule(out[i]));
           i += j; // skip empties
        }
     }

     // checking for real stoichiometric coeffs
     bool real(false);
     for(unsigned int i = 0; i < _reactants.size(); i++)
     {
        if(this->after_coma_digits((_reactants[i].second)))
        {
           real = true;
           break;
        }
     }
     if(!real)
     {
       for(unsigned int i = 0; i < _products.size(); i++)
       {
         if(this->after_coma_digits(_products[i].second))
         {
           real = true;
           break;
         }
       }
    }

    if(real)this->rescale_stoichiometry();
  }

  template <typename NumericType>
  inline
  std::pair<std::string,NumericType> ChemKinParser<NumericType>::parse_molecule(const std::string & molecule)
  {
     // as long as we have a number ({[0-9],.}), it is the stoichiometric coefficient
     unsigned int pos(0);
     while(this->is_real_number(molecule[pos]))pos++;

     NumericType stoi = (pos == 0)?1.:std::atof(molecule.substr(0,pos+1).c_str());

     return std::make_pair(molecule.substr(pos,std::string::npos),stoi);
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::is_real_number(const char & c) const
  {
     return (c == '0' || c == '1' || c == '2' ||
             c == '3' || c == '4' || c == '5' ||
             c == '6' || c == '7' || c == '8' ||
             c == '9' || c == '.');
             
  }

  template <typename NumericType>
  inline
  void ChemKinParser<NumericType>::parse_coefficients_line(const std::string &line)
  {
      std::vector<std::string> out;
      int nk = SplitString(line,_spec.parser(),out,false);
      if(nk < 2)antioch_parsing_error("ChemKin parser: can't parse this line:\n" + line);
      // can be LOW, TROE or coefficients
        //accounts for blank spaces
      if(out.front().find(_map.at(ParsingKey::TROE_FALLOFF)) != std::string::npos) //TROE, alpha, T***, T*, T**
      {
        antioch_assert_greater_equal(out.size(),2);
        std::vector<std::string> troe_par;
        int npar = SplitString(out[1]," ",troe_par,false);
        if(npar < 3)antioch_parsing_error("ChemKin parser: Troe parameters error while reading:\n" + line);

         _Troe_alpha = std::atof(troe_par[0].c_str());
         _Troe_T3    = std::atof(troe_par[1].c_str());
         _Troe_T1    = std::atof(troe_par[2].c_str());
         _Troe_T2    = (npar == 4)?std::atof(troe_par[3].c_str()):-1.L;
      }
      else if(out.front().find(_map.at(ParsingKey::FALLOFF_LOW_NAME)) != std::string::npos) // k0
      {
        antioch_assert_greater_equal(out.size(),2);
        std::vector<std::string> k0_par;
        int npar = SplitString(out[1]," ",k0_par,false);
        if(npar != 3)antioch_parsing_error("ChemKin parser: Falloff k0 parameters error while reading:\n" + line);

        _A.insert(_A.begin(),std::atof(k0_par[0].c_str()));
        _b.insert(_b.begin(),std::atof(k0_par[1].c_str()));
        _Ea.insert(_Ea.begin(),std::atof(k0_par[2].c_str()));
         
      }else // efficiencies
      {
        for(unsigned int i = 0; i < out.size(); i++)
        {
            // Trim from the left
            out[i].erase(out[i].begin(), std::find_if(out[i].begin(), out[i].end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
            // Trim from the right
            out[i].erase(std::find_if(out[i].rbegin(), out[i].rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), out[i].end());

            if(out[i].empty())
            {
              out.erase(out.begin() + i);
              i--;
            }
        }
        if(out.size()%2 != 0)antioch_parsing_error("ChemKin parser: efficiencies parsing failed:\n" + line);
        for(unsigned int c = 0; c < out.size(); c += 2)
        {
           _efficiencies.push_back(std::make_pair( out[c], std::atof(out[c+1].c_str()))); // mol / coef
        }
      }
      
  }

  template <typename NumericType>
  inline
  void ChemKinParser<NumericType>::rescale_stoichiometry()
  {
// we suppose rational number r = p / q with p and q integers, 
// we look for q kinda stupidly (if this is more than 2 or 3, it
// is really ridiculous)
     std::vector<unsigned int> not_int_factor;
     for(unsigned int i = 0; i < _reactants.size(); i++)
     {
       if(this->after_coma_digits(_reactants[i].second))
       {
         unsigned int fac = this->factor_to_int(_reactants[i].second);
         unsigned int i(0);
         for(i = 0; i < not_int_factor.size(); i++)
         {
           if(fac == not_int_factor[i])break;
         }
         if(i < not_int_factor.size() - 1)not_int_factor.push_back(fac);
       }
     }
     for(unsigned int i = 0; i < _products.size(); i++)
     {
       if(this->after_coma_digits(_products[i].second))
       {
         unsigned int fac = this->factor_to_int(_products[i].second);
         unsigned int i(0);
         for(i = 0; i < not_int_factor.size(); i++)
         {
           if(fac == not_int_factor[i])break;
         }
         if(i < not_int_factor.size() - 1)not_int_factor.push_back(fac);
       }
     }

// global factor will be brutal multiplication
// don't want to decompose into prime numbers
     
     unsigned int factor(1);
     for(unsigned int i = 0; i < not_int_factor.size(); i++)
     {
         factor *= not_int_factor[i];
     }

// rescale everyone
     for(unsigned int i = 0; i < _reactants.size(); i++)
     {
        _reactants[i].second *= (NumericType)factor;
     }
     for(unsigned int i = 0; i < _products.size(); i++)
     {
        _products[i].second *= (NumericType)factor;
     }
     
 }

 template <typename NumericType>
 inline
 unsigned int ChemKinParser<NumericType>::factor_to_int(NumericType number) const
 {
    // now find p and q associated with this number
    unsigned int down(2);
    const unsigned int limit(150);
    while(this->after_coma_digits(number * (NumericType) down))
    {
      down++;
      if(down > limit)
      {
         std::stringstream os;
         os << "real is " << number << " and multiplicative factor limit is " << limit;
         antioch_parsing_error("ChemKin parser: could not find integer from real\n:" + os.str());
      }
    }

    return down;
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::after_coma_digits(NumericType number) const
  {
    // smallest number tolerated 1e-3, it is already ridiculous
    const NumericType eps(1e-3);
    return !((number - std::floor(number)) > eps);
  }

}//end namespace Antioch

#endif
