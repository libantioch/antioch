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
#ifndef ANTIOCH_XML_PARSER_H
#define ANTIOCH_XML_PARSER_H

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/string_utils.h"
#include "antioch/parser_base.h"
#include "antioch/parsing_enum.h"

//C++
#include <string>
#include <vector>
#include <map>

namespace tinyxml2
{
  class XMLDocument;
  class XMLElement;
}

namespace Antioch{


  template <typename CoeffType>
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

  /*!\class XMLParser

     Nothing is stored, this parser is based on the tinyxml2
     implementation. Please note that no other file should include
     the `tinyxml2_imp.h' header.

     The defaults units are based and derived on Cantera:
       -   pre-exponential parameters in (m3/kmol)^(m-1)/s
       -   activation energy in cal/mol,
       -   power parameter without unit
       -   cross-section typically in cm2/nm,
       -   lambda typically in nm,
   */
  template <typename NumericType = double>
  class XMLParser: public ParserBase<NumericType>
  {
        public:
          XMLParser(const std::string &filename, bool verbose = true);
          ~XMLParser();

          void change_file(const std::string & filename);

//// first local pointers
         /*! Read header of file, go to interesting part*/
         bool initialize();

/// species
        //! reads the species set
        const std::vector<std::string> species_list() ;

        //! reads the thermo, NASA generalist, no templates for virtual
        void read_thermodynamic_data(NASAThermoMixture<NumericType, NASA7CurveFit<NumericType> >& thermo)
                {this->read_thermodynamic_data_root(thermo);}

        //! reads the thermo, NASA generalist, no templates for virtual
        void read_thermodynamic_data(NASAThermoMixture<NumericType, NASA9CurveFit<NumericType> >& thermo)
                {this->read_thermodynamic_data_root(thermo);}

        //! reads the thermo, NASA generalist, no templates for virtual
        void read_thermodynamic_data(NASAThermoMixture<NumericType, CEACurveFit<NumericType> >& thermo)
                {this->read_thermodynamic_data_root(thermo);}

        //! reads the thermo, CEA deprecated
        void read_thermodynamic_data(CEAThermodynamics<NumericType >& thermo)
                {this->read_thermodynamic_data_root(thermo);}

/// reaction

         /*! go to next reaction*/
         bool reaction();

         /*! go to next rate constant*/
         bool rate_constant(const std::string & kinetics_model);

         /*! return true if there's a Troe block*/
         bool Troe() const;

         /*! return reaction id, 0 if not provided*/
         const std::string reaction_id() const;

         /*! return reaction equation */
         const std::string reaction_equation() const;

         /*! return reaction chemical process*/
         const std::string reaction_chemical_process() const;

         /*! return reaction kinetics model*/
         const std::string reaction_kinetics_model(const std::vector<std::string> &kinetics_models) const;

         /*! return reversible state*/
         bool reaction_reversible() const;

         /*! return pairs of reactants and stoichiometric coefficients*/
         bool reactants_pairs(std::vector<std::pair<std::string,int> > & reactants_pair) const;

         /*! return pairs of products and stoichiometric coefficients*/
         bool products_pairs(std::vector<std::pair<std::string,int> > & products_pair) const;

         /*! return a map between reactants' name and found partial orders */
         const std::map<std::string,NumericType> reactants_orders() const;

         /*! return a map between products' name and found partial orders */
         const std::map<std::string,NumericType> products_orders() const;

         /*! return true if "name" attribute is found with value "k0"*/
         bool is_k0(unsigned int nrc, const std::string & kin_model) const;

         /*! return index of k0 (0 or 1)*/
         unsigned int where_is_k0(const std::string & kin_model) const;

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
         template <typename ThermoType>
         void read_thermodynamic_data_root(ThermoType & thermo);

         /*! return pairs of molecules and stoichiometric coefficients*/
         template <typename PairedType>
         bool molecules_pairs(tinyxml2::XMLElement * molecules, std::vector<std::pair<std::string,PairedType> > & products_pair) const;

         /*! return a parameter's value*/
         bool get_parameter(const tinyxml2::XMLElement * ptr, const std::string & par, NumericType & par_value, std::string & par_unit) const;

         /*! return a parameter's values*/
         bool get_parameter(const tinyxml2::XMLElement * ptr, const std::string & par, std::vector<NumericType> & numpar, std::string & par_unit) const;

         /*! return the unit of current pointer*/
         const std::string unit(tinyxml2::XMLElement * parameter) const;


          /*! Never use default constructor*/
          XMLParser();
          tinyxml2::XMLDocument * _doc;

//
          tinyxml2::XMLElement * _species_block;
          tinyxml2::XMLElement * _thermo_block;
          tinyxml2::XMLElement * _reaction_block;
          tinyxml2::XMLElement * _reaction;

          tinyxml2::XMLElement * _rate_constant;
          tinyxml2::XMLElement * _Troe;

          std::map<ParsingKey,std::string> _map;
          std::map<ParsingKey,std::string> _default_unit;

  };

}//end namespace Antioch

#endif
