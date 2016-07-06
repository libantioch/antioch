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

#ifndef ANTIOCH_PARSER_BASE_H
#define ANTIOCH_PARSER_BASE_H

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/parsing_enum.h"
#include "antioch/input_utils.h"

//C++
#include <vector>
#include <string>
#include <sstream>
#include <map>

namespace Antioch
{

  // Forward declarations
  template <class NumericType>
  class ChemicalMixture;

  template <class NumericType>
  class TransportMixture;

  template <typename NumericType, typename CurveFit>
  class NASAEvaluator;

  template <typename NumericType, typename CurveFit>
  class NASAThermoMixture;

  template <typename NumericType>
  class NASA7CurveFit;

  template <typename NumericType>
  class NASA9CurveFit;

  // backward compatibility
  template <typename NumericType>
  class CEACurveFit;

  template <typename NumericType>
  class CEAEvaluator;

// micro
  template <typename NumericType>
  class StatMechThermodynamics;

  template <typename Macro, typename NumericType>
  class IdealGasMicroThermo;

  /*!\class ParserBase

      A parser is an instance related to a file. The parser
      corresponds to a file type (e.g. XML or ChemKin). The
      file HAS to be given in the constructor as the parser
      is indissociable from the file. Set to `true' by default,
      a verbose switch is also available.

      We define here the rule of parsing for
      all parsers. Differences/specificities are described in the
      corresponding files.

      The different things to parse are:
        - the species, name and core characteristics:
                - list of names (std::string) `const std::vector<std::string> species_list() const;'
                - mandatory data (molar mass, heat of formation at 0 K,
                                  number of translational/rotational DOFs,
                                  charge number) `void read_chemical_species(ChemicalMixture<NumericType> & chem_mixture);'
                - vibrational data `void read_vibrational_data(ChemicalMixture<NumericType> & chem_mixture);'
                - electronic data `void read_electronic_data(ChemicalMixture<NumericType> & chem_mixture);'
        - the kinetics:
                - it requires a boolean to ensure there is a reaction to parse `bool reaction();'
                - ...
   */

  template <typename NumericType>
  class ParserBase
  {
  public:
    ParserBase(const std::string & type,
               const std::string & file,
               bool verbose = true,
               const std::string & comments = "");

    virtual ~ParserBase(){};

        // initialize kinetics, mandatory
        virtual bool initialize() = 0;

        // to reinitialize, mandatory
        virtual void change_file(const std::string & filename) = 0;

/// species
        //! reads the species set
        virtual const std::vector<std::string> species_list() {antioch_not_implemented_msg(_not_implemented); return std::vector<std::string>();}

        //! reads the mandatory data, not valid in xml && chemkin
        virtual void read_chemical_species(ChemicalMixture<NumericType> & /*chem_mixture*/)  {antioch_not_implemented_msg(_not_implemented);}

        //! reads the vibrational data, not valid in xml && chemkin
        virtual void read_vibrational_data(ChemicalMixture<NumericType> & /*chem_mixture*/)  {antioch_not_implemented_msg(_not_implemented);}

        //! reads the electronic data, not valid in xml && chemkin
        virtual void read_electronic_data(ChemicalMixture<NumericType> & /*chem_mixture*/)  {antioch_not_implemented_msg(_not_implemented);}

// transport

        //! reads the transport data, not valid in xml && chemkin
        virtual void read_transport_data(TransportMixture<NumericType> & /*transport_mixture*/)  {antioch_not_implemented_msg(_not_implemented);}

/// thermo

        //! reads the thermo, NASA generalist, no templates for virtual
        virtual void read_thermodynamic_data(NASAThermoMixture<NumericType, NASA7CurveFit<NumericType> >& /*thermo*/)  {antioch_not_implemented_msg(_not_implemented);}

        //! reads the thermo, NASA generalist, no templates for virtual
        virtual void read_thermodynamic_data(NASAThermoMixture<NumericType, NASA9CurveFit<NumericType> >& /*thermo*/)  {antioch_not_implemented_msg(_not_implemented);}

        //! reads the thermo, NASA generalist, no templates for virtual
        virtual void read_thermodynamic_data(NASAThermoMixture<NumericType, CEACurveFit<NumericType> >& /*thermo*/)  {antioch_not_implemented_msg(_not_implemented);}


/// reaction

// non const

         /*! read & store current reaction and go to next reaction*/
         virtual bool reaction() {antioch_not_implemented_msg(_not_implemented); return false;}

         /*! go to next rate constant*/
         virtual bool rate_constant(const std::string & /*kinetics_model*/) {antioch_not_implemented_msg(_not_implemented); return false;}

// const

         /*! \return true if there's a Troe block*/
         virtual bool Troe() const {antioch_not_implemented_msg(_not_implemented); return false;}

         /*! \return reaction id, 0 if not provided*/
         virtual const std::string reaction_id() const  {antioch_not_implemented_msg(_not_implemented); return std::string();}

         /*! \return reaction equation */
         virtual const std::string reaction_equation() const {antioch_not_implemented_msg(_not_implemented); return std::string();}

         /*! \return reaction chemical process*/
         virtual const std::string reaction_chemical_process() const {antioch_not_implemented_msg(_not_implemented); return std::string();}

         /*! \return reaction kinetics model*/
         virtual const std::string reaction_kinetics_model(const std::vector<std::string> & /*kinetics_models*/) const  {antioch_not_implemented_msg(_not_implemented); return std::string();}

         /*! \return reversible state*/
         virtual bool reaction_reversible() const {antioch_not_implemented_msg(_not_implemented); return false;}

         /*! \return pairs of reactants and stoichiometric coefficients*/
         virtual bool reactants_pairs(std::vector<std::pair<std::string,int> > & /*reactants_pair*/) const  {antioch_not_implemented_msg(_not_implemented); return false;}

         /*! \return pairs of products and stoichiometric coefficients*/
         virtual bool products_pairs(std::vector<std::pair<std::string,int> > & /*products_pair*/) const {antioch_not_implemented_msg(_not_implemented); return false;}

         /*! return a map between reactants' name and found partial orders */
         virtual const std::map<std::string,NumericType> reactants_orders() const {antioch_not_implemented_msg(_not_implemented); return std::map<std::string,NumericType>();}

         /*! return a map between products' name and found partial orders */
         virtual const std::map<std::string,NumericType> products_orders() const {antioch_not_implemented_msg(_not_implemented); return std::map<std::string,NumericType>();}

         /*! \return true if "name" attribute is found with value "k0"*/
         virtual bool is_k0(unsigned int /*nrc*/, const std::string & /*kin_model*/) const {antioch_not_implemented_msg(_not_implemented); return false;}

         /*! \return index of k0 (0 or 1)*/
         virtual unsigned int where_is_k0(const std::string & /*kin_model*/) const {antioch_not_implemented_msg(_not_implemented); return -1;}

         /*! \return true if pre exponentiel coefficient*/
         virtual bool rate_constant_preexponential_parameter(NumericType & /*A*/, std::string & /*A_unit*/, std::string & /*def_unit*/) const {antioch_not_implemented_msg(_not_implemented); return false;}

         /*! \return true if beta coefficient*/
         virtual bool rate_constant_power_parameter(NumericType & /*b*/, std::string & /*b_unit*/, std::string & /*def_unit*/) const {antioch_not_implemented_msg(_not_implemented); return false;}

         /*! \return true if activation energie*/
         virtual bool rate_constant_activation_energy_parameter(NumericType & /*Ea*/, std::string & /*Ea_unit*/, std::string & /*def_unit*/) const {antioch_not_implemented_msg(_not_implemented); return false;}

         /*! \return true if D coefficient*/
         virtual bool rate_constant_Berthelot_coefficient_parameter(NumericType & /*D*/, std::string & /*D_unit*/, std::string & /*def_unit*/) const {antioch_not_implemented_msg(_not_implemented); return false;}

         /*! \return true if Tref*/
         virtual bool rate_constant_Tref_parameter( NumericType & /*Tref*/, std::string & /*Tref_unit*/, std::string & /*def_unit*/) const {antioch_not_implemented_msg(_not_implemented); return false;}

         /*! \return true if lambda*/
         virtual bool rate_constant_lambda_parameter(std::vector<NumericType> & /*lambda*/, std::string & /*lambda_unit*/, std::string & /*def_unit*/) const {antioch_not_implemented_msg(_not_implemented); return false;}

         /*! \return true if sigma*/
         virtual bool rate_constant_cross_section_parameter(std::vector<NumericType> & /*sigma*/,  std::string & /*sigma_unit*/, std::string & /*def_unit*/) const {antioch_not_implemented_msg(_not_implemented); return false;}

         /*! \return true if a Kooij is called Arrhenuis*/
         virtual bool verify_Kooij_in_place_of_Arrhenius() const {antioch_not_implemented_msg(_not_implemented); return false;}

         /*! \return true if efficiencies are found*/
         virtual bool efficiencies(std::vector<std::pair<std::string,NumericType> > & /*par_values*/) const {antioch_not_implemented_msg(_not_implemented); return false;}

         /*! \return true is alpha*/
         virtual bool Troe_alpha_parameter(NumericType & /*alpha*/, std::string & /*alpha_unit*/, std::string & /*def_unit*/) const {antioch_not_implemented_msg(_not_implemented); return false;}

         /*! \return true is alpha*/
         virtual bool Troe_T1_parameter(NumericType & /*T1*/, std::string & /*T1_unit*/, std::string & /*def_unit*/) const {antioch_not_implemented_msg(_not_implemented); return false;}

         /*! \return true is alpha*/
         virtual bool Troe_T2_parameter(NumericType & /*T2*/, std::string & /*T2_unit*/, std::string & /*def_unit*/) const {antioch_not_implemented_msg(_not_implemented); return false;}

         /*! \return true is alpha*/
         virtual bool Troe_T3_parameter(NumericType & /*T3*/, std::string & /*T3_unit*/, std::string & /*def_unit*/) const {antioch_not_implemented_msg(_not_implemented); return false;}

         /*! \return name of file*/
        const std::string file() const {return _file;}

         /*! \return type of parser*/
        const std::string type() const {return _type;}

         /*! \return verbosity*/
         bool verbose() const {return _verbose;}

        /*! \return enum*/
        ParsingType enum_type() const;


     protected:

        void skip_comments(std::istream & doc);

        std::string _type;
        std::string _file;
        bool        _verbose;
        std::string _comments;

        std::string _not_implemented;

     private:
        ParserBase();
  };

} // end namespace Antioch


#endif //ANTIOCH_PARSER_BASE_H
