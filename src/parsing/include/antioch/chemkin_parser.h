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
#ifndef ANTIOCH_CHEMKIN_PARSER_H
#define ANTIOCH_CHEMKIN_PARSER_H

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/string_utils.h"
#include "antioch/parser_base.h"
#include "antioch/parsing_enum.h"
#include "antioch/units.h"
#include "antioch/chemkin_definitions.h"
#include "antioch/nasa_mixture.h" // because not templated, therefore should be entirely known
#include "antioch/cea_curve_fit.h" // because not templated, therefore should be entirely known
#include "antioch/nasa_curve_fit.h" // because not templated, therefore should be entirely known

//ChemKin

//C++
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

namespace Antioch{


  template <typename NumericType>
  class ChemicalMixture;

  template <typename NumericType, typename CurveType>
  class NASAThermoMixture;

  template <typename NumericType>
  class NASA7CurveFit;

  template <typename NumericType>
  class NASA9CurveFit;

  // backward compatibility
  template <typename NumericType>
  class CEACurveFit;

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
  class ChemKinParser: public ParserBase<NumericType>
  {
        public:
          ChemKinParser(const std::string &filename, bool verbose = true);
          ~ChemKinParser();

         void change_file(const std::string & filename);

////////////// species

         /*! read SPECIES block*/
         const std::vector<std::string> species_list();

        // reads the mandatory data, not valid yet in ChemKin
//        void read_chemical_species(ChemicalMixture<NumericType> & chem_mixture);

        // reads the vibrational data, not valid yet in ChemKin
//        void read_vibrational_data(ChemicalMixture<NumericType> & chem_mixture);

        // reads the electronic data, not valid yet in ChemKin
//        void read_electronic_data(ChemicalMixture<NumericType> & chem_mixture);

////////////////// thermo

//global overload
// it seems that they're all linked in
// a way so they shadow themselves
// => we need to implement all or nothing

        //! reads the thermo, NASA generalist, no templates for virtual
        void read_thermodynamic_data(NASAThermoMixture<NumericType, NASA7CurveFit<NumericType> >& thermo)
                {this->read_thermodynamic_data_root(thermo);}

        //! reads the thermo, NASA generalist, no templates for virtual
        void read_thermodynamic_data(NASAThermoMixture<NumericType, NASA9CurveFit<NumericType> >& /*thermo*/)
                {antioch_not_implemented_msg(ParserBase<NumericType>::_not_implemented);}

        //! reads the thermo, NASA generalist, no templates for virtual
        void read_thermodynamic_data(NASAThermoMixture<NumericType, CEACurveFit<NumericType> >& /*thermo*/)
                {antioch_not_implemented_msg(ParserBase<NumericType>::_not_implemented);}

        //! reads the thermo, CEA deprecated
        void read_thermodynamic_data(CEAThermodynamics<NumericType >& /*thermo*/)
                {antioch_not_implemented_msg(ParserBase<NumericType>::_not_implemented);}

///////////////// kinetics

         /*! Read header of file, go to interesting part*/
         bool initialize();

         /*! go to next reaction*/
         bool reaction();

         /*! go to next rate constant*/
         bool rate_constant(const std::string & /* kinetics_model */);

         /*! return true if there's a Troe block*/
         bool Troe() const;

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

          //! reads the thermo, NASA generalist
          template <typename CurveType>
          void read_thermodynamic_data_root(NASAThermoMixture<NumericType, CurveType >& thermo);

          /*! Convenient method */
          void parse_a_line(const std::string & line);

          /*! Convenient method */
          void parse_equation_coef(const std::string & line);

          /*! Convenient method */
          void parse_reversible_parameters(const std::string & line);

          /*! Convenient method */
          void parse_equation(std::string &equation);

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

          /*! verify if line is a new reaction*/
          bool next_reaction(const std::string & line);

          /*! finding next line that might be a reaction */
          bool next_meaningful_line(std::string & line);

          /*! Never use default constructor*/
          ChemKinParser();
          std::ifstream                    _doc;


          bool                             _reversible;
          unsigned int                     _nrates; // total number of rates
          unsigned int                     _crates; // current place

          unsigned int                     _pow_unit; // for A unit

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

          std::string _cached_line;
          bool        _duplicate_process;
          bool        _next_is_reverse;

          const ChemKinDefinitions _spec;

  };

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::Troe() const
  {
      return (_chemical_process.find("TroeFalloff") != std::string::npos);
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
  bool ChemKinParser<NumericType>::is_k0(unsigned int nrc, const std::string & /*kin_model*/) const
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
  bool ChemKinParser<NumericType>::rate_constant_Berthelot_coefficient_parameter(NumericType & /*D*/, std::string & /*D_unit*/, std::string & /*def_unit*/) const
  {
      return false; //not a chemkin model
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::rate_constant_lambda_parameter(std::vector<NumericType> & /*lambda*/, std::string & /*lambda_unit*/, std::string & /*def_unit*/) const
  {
      return false; //not a supported chemkin model yet
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::rate_constant_cross_section_parameter(std::vector<NumericType> & /*sigma*/, std::string & /*sigma_unit*/, std::string & /*def_unit*/) const
  {
      return false; //not a supported chemkin model yet
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
  template <typename CurveType>
  inline
  void ChemKinParser<NumericType>::read_thermodynamic_data_root(NASAThermoMixture<NumericType, CurveType >& thermo)
  {

// finding thermo
    std::string line;
    ascii_getline(_doc,line);

    while(line.find(_map.at(ParsingKey::THERMO)) == std::string::npos)
    {
       if(!ascii_getline(_doc,line) || _doc.eof())break;
    }

    if(!_doc.good())
    {
        std::cerr << "Thermodynamics description not found" << std::endl;
        antioch_error();
    }

    const ChemicalMixture<NumericType>& chem_mixture = thermo.chemical_mixture();

    std::string name;
    std::vector<NumericType> coeffs;
    std::vector<NumericType> temps(3,0.);

    this->skip_comments(_doc);

// chemkin classic
// only two intervals
// \todo: the parser should allow custom
// intervals definition
    while (_doc.good())
      {
        std::stringstream tmp;
        this->skip_comments(_doc); // comments in middle

        if(!ascii_getline(_doc,line))break;

        if(line.find(_spec.end_tag()) != std::string::npos)break; //END

        // question is: will we have temps or NASA coeffs
        // line is 80 character long for coeffs, so if not, it's temps
        if(line.size() != 80)
        {
           std::vector<std::string> temp_tmp;
           SplitString(line," ",temp_tmp,false);
           if(temp_tmp.size() == 0)
           {
              std::cerr << "This line is not understandable by the chemkin parser:\n" << line << std::endl;
              antioch_error();
           }
           temps.resize(temp_tmp.size(),0.);
           for(unsigned int t = 0; t < temp_tmp.size(); t++)
           {
               temps[t] = std::atof(temp_tmp[t].c_str());
           }
           continue;
        }

        tmp << line.substr(0,18); //name

        // this is ChemKin doc, 1st: temps E10.0
        tmp << line.substr(45,10);         // low
        tmp << " " << line.substr(55,10);  // high
        tmp << " " << line.substr(65,10);  // inter

         // get rid of last character
        for(unsigned int n = 0; n < 3; n++)
        {
          if(!ascii_getline(_doc,line))// we have a problem
          {
             std::cerr << "NASA input file error" << std::endl;
             antioch_error();
          }
                //2nd: coeffs E15.8
          tmp << line.substr(0,15)  << " "
              << line.substr(15,15) << " "
              << line.substr(30,15) << " "
              << line.substr(45,15) << " "
              << line.substr(60,15) << " ";
        }

        coeffs.clear();
        tmp >> name     // Species Name
            >> temps[0]
            >> temps[2]
            >> temps[1];
        for(unsigned int i = 0; i < 14; i++)
        {
           NumericType a;
           tmp >> a;
           coeffs.push_back(a);
        }

        // If we are still good, we have a valid set of thermodynamic
        // data for this species. Otherwise, we read past end-of-file
        // in the section above
        if (_doc.good())
          {
            // Check if this is a species we want.
            if( chem_mixture.species_name_map().find(name) !=
                chem_mixture.species_name_map().end() )
              {
                thermo.add_curve_fit(name, coeffs, temps);
              }
          }
      } // end while
  }

}//end namespace Antioch

#endif
