//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_PARSER_BASE_H
#define ANTIOCH_PARSER_BASE_H

//Antioch
#include "antioch/parsing_enum.h"

//C++
#include <vector>
#include <string>

namespace Antioch
{

  // Forward declarations
  template <class NumericType>
  class ChemicalMixture;

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
        ParserBase(const std::string & type, const std::string & file, bool verbose = true);
        virtual ~ParserBase();

        // to reinitialize, mandatory
        virtual void change_file(const std::string & filename) = 0;

/// species
        //! reads the species set
        virtual const std::vector<std::string> species_list() const {not_implemented(); return std::vector<std::string>();}

        //! reads the mandatory data, not valid in xml
        virtual void read_chemical_species(ChemicalMixture<NumericType> & /*chem_mixture*/)  {not_implemented();}

        //! reads the vibrational data, not valid in xml
        virtual void read_vibrational_data(ChemicalMixture<NumericType> & /*chem_mixture*/)  {not_implemented();}

        //! reads the electronic data, not valid in xml
        virtual void read_electronic_data(ChemicalMixture<NumericType> & /*chem_mixture*/)  {not_implemented();}

/// reaction

         /*! go to next reaction*/
         virtual bool reaction() {not_implemented(); return false;}

         /*! go to next rate constant*/
         virtual bool rate_constant(const std::string & /*kinetics_model*/) {not_implemented(); return false;}

         /*! \return true if there's a Troe block*/
         virtual bool Troe()  {not_implemented(); return false;}

         /*! \return reaction id, 0 if not provided*/
         virtual const std::string reaction_id() const  {not_implemented(); return std::string();}

         /*! \return reaction equation */
         virtual const std::string reaction_equation() const {not_implemented(); return std::string();}

         /*! \return reaction chemical process*/
         virtual const std::string reaction_chemical_process() const {not_implemented(); return std::string();}

         /*! \return reaction kinetics model*/
         virtual const std::string reaction_kinetics_model(const std::vector<std::string> & /*kinetics_models*/) const  {not_implemented(); return std::string();}

         /*! \return reversible state*/
         virtual bool reaction_reversible() const {not_implemented(); return false;}

         /*! \return pairs of reactants and stoichiometric coefficients*/
         virtual bool reactants_pairs(std::vector<std::pair<std::string,int> > & /*reactants_pair*/) const  {not_implemented(); return false;}

         /*! \return pairs of products and stoichiometric coefficients*/
         virtual bool products_pairs(std::vector<std::pair<std::string,int> > & /*products_pair*/) const {not_implemented(); return false;}

         /*! \return true if "name" attribute is found with value "k0"*/
         virtual bool is_k0(unsigned int nrc, const std::string & /*kin_model*/) const {not_implemented(); return false;}

         /*! \return index of k0 (0 or 1)*/
         virtual unsigned int where_is_k0(const std::string & /*kin_model*/) const {not_implemented(); return -1;}

         /*! \return true if pre exponentiel coefficient*/
         virtual bool rate_constant_preexponential_parameter(NumericType & /*A*/, std::string & /*A_unit*/, std::string & /*def_unit*/) const {not_implemented(); return false;}

         /*! \return true if beta coefficient*/
         virtual bool rate_constant_power_parameter(NumericType & /*b*/, std::string & /*b_unit*/, std::string & /*def_unit*/) const {not_implemented(); return false;}

         /*! \return true if activation energie*/
         virtual bool rate_constant_activation_energy_parameter(NumericType & /*Ea*/, std::string & /*Ea_unit*/, std::string & /*def_unit*/) const {not_implemented(); return false;}

         /*! \return true if D coefficient*/
         virtual bool rate_constant_Berthelot_coefficient_parameter(NumericType & /*D*/, std::string & /*D_unit*/, std::string & /*def_unit*/) const {not_implemented(); return false;}

         /*! \return true if Tref*/
         virtual bool rate_constant_Tref_parameter( NumericType & /*Tref*/, std::string & /*Tref_unit*/, std::string & /*def_unit*/) const {not_implemented(); return false;}

         /*! \return true if lambda*/
         virtual bool rate_constant_lambda_parameter(std::vector<NumericType> & /*lambda*/, std::string & /*lambda_unit*/, std::string & /*def_unit*/) const {not_implemented(); return false;}

         /*! \return true if sigma*/
         virtual bool rate_constant_cross_section_parameter(std::vector<NumericType> & /*sigma*/,  std::string & /*sigma_unit*/, std::string & /*def_unit*/) const {not_implemented(); return false;}

         /*! \return true if a Kooij is called Arrhenuis*/
         virtual bool verify_Kooij_in_place_of_Arrhenius() const {not_implemented(); return false;}

         /*! \return true if efficiencies are found*/
         virtual bool efficiencies(std::vector<std::pair<std::string,NumericType> > & /*par_values*/) const {not_implemented(); return false;}

         /*! \return true is alpha*/
         virtual bool Troe_alpha_parameter(NumericType & /*alpha*/, std::string & /*alpha_unit*/, std::string & /*def_unit*/) const {not_implemented(); return false;}

         /*! \return true is alpha*/
         virtual bool Troe_T1_parameter(NumericType & /*T1*/, std::string & /*T1_unit*/, std::string & /*def_unit*/) const {not_implemented(); return false;}

         /*! \return true is alpha*/
         virtual bool Troe_T2_parameter(NumericType & /*T2*/, std::string & /*T2_unit*/, std::string & /*def_unit*/) const {not_implemented(); return false;}

         /*! \return true is alpha*/
         virtual bool Troe_T3_parameter(NumericType & /*T3*/, std::string & /*T3_unit*/, std::string & /*def_unit*/) const {not_implemented(); return false;}

         /*! \return name of file*/
        const std::string file() const {return _file;}

         /*! \return type of parser*/
        const std::string type() const {return _type;}

         /*! \return verbosity*/
        const bool     verbose() const {return _verbose;}

        /*! \return enum*/
        ParsingType enum_type() const;
        

     protected:

        std::string _type;
        std::string _file;
        bool        _verbose;

        void not_implemented() const;

     private:
        ParserBase();
  };

  template <typename NumericType>
  ParserBase<NumericType>::ParserBase(const std::string & type, const std::string & file, bool verbose):
        _type(type),
        _file(file),
        _verbose(verbose)
  {
      return;
  }

  template <typename NumericType>
  ParserBase<NumericType>::~ParserBase()
  {
      return;
  }

  template <typename NumericType>
  ParsingType ParserBase<NumericType>::enum_type() const
  {
      ParsingType PType(ASCII);
      if(_type == "ascii")
      {
          PType = ASCII;
      }else if(_type == "ChemKin")
      {
          PType = CHEMKIN;
      }else if(_type == "XML")
      {
          PType = XML;
      }else
      {
           antioch_parsing_error(std::string("unknown parser type!!! " + _type));
      }
      return PType;
  }

  template <typename NumericType>
  void ParserBase<NumericType>::not_implemented() const
  {
      std::cerr << "\n*********************************************************\n"
                << "This method is not available with a " << _type << " parser.\n"
                << "Parsing file " << _file << ".\n"
                << "No format has been defined yet.  Maybe contribute?\n"
                << "https://github.com/libantioch/antioch\n" 
                << "\n*********************************************************\n"
                << std::endl;
  }


} // end namespace Antioch


#endif //ANTIOCH_PARSER_BASE_H
