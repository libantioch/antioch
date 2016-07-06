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

#ifndef ANTIOCH_ASCII_PARSER_H
#define ANTIOCH_ASCII_PARSER_H

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/parser_base.h"
#include "antioch/parsing_enum.h"
#include "antioch/input_utils.h"
#include "antioch/string_utils.h"
#include "antioch/units.h"

// C++
#include <fstream>
#include <string>
#include <vector>
#include <map>

namespace Antioch
{

   // backward compatibility
  typedef unsigned int Species;

  // Forward declarations
  template <class NumericType>
  class ChemicalMixture;

  template <class NumericType>
  class TransportMixture;

// macro
  template <typename NumericType, typename CurveFit>
  class NASAThermoMixture;

  template <typename NumericType, typename CurveFit>
  class NASAEvaluator;

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

  template <typename NumericType>
  class ASCIIParser: public ParserBase<NumericType>
  {
      public:
        ASCIIParser(const std::string& file, bool verbose = true);
        ~ASCIIParser();

        void change_file(const std::string & filename);

        bool initialize() {return false;}

        //! set the indexes of to-be-ignored columns
        void set_ignored_columns(const std::vector<unsigned int> & ignored);

/////////////////
// species
////////////////
        //! read species list
        const std::vector<std::string>  species_list() ;

        //! read the mandatory data
        void read_chemical_species(ChemicalMixture<NumericType> & chem_mixture);

        //! read the vibrational data
        void read_vibrational_data(ChemicalMixture<NumericType>& chem_mixture);

        //! read the electronic data
        void read_electronic_data(ChemicalMixture<NumericType>& chem_mixture);

        //! reads the transport data, not valid in xml && chemkin
        void read_transport_data(TransportMixture<NumericType> & transport_mixture) {this->read_transport_data_root(transport_mixture);}


///////////////
// thermo
///////////////

//global overload
    //! reads the thermo, NASA generalist, no templates for virtual
    void read_thermodynamic_data(NASAThermoMixture<NumericType, NASA7CurveFit<NumericType> >& /*thermo*/)
    {antioch_error_msg("ERROR: ASCIIParsing only supports parsing for CEACurveFit!");}

    //! reads the thermo, NASA generalist, no templates for virtual
    void read_thermodynamic_data(NASAThermoMixture<NumericType, NASA9CurveFit<NumericType> >& /*thermo*/)
    {antioch_error_msg("ERROR: ASCIIParsing only supports parsing for CEACurveFit!");}

        //! reads the thermo, NASA generalist, no templates for virtual
    void read_thermodynamic_data(NASAThermoMixture<NumericType, CEACurveFit<NumericType> >& thermo)
    {this->read_thermodynamic_data_root(thermo);}

     private:

// templated thermo version
        //! read the thermodynamic data
        template <typename CurveType>
        void read_thermodynamic_data_root(NASAThermoMixture<NumericType, CurveType >& thermo);

// templated transport version
        //! read the thermodynamic data
        template <typename Mixture>
        void read_transport_data_root(Mixture & transport);

        //! find the index of the wanted data
        void find_first(unsigned int & index,unsigned int n_data) const;

    //! Warn about species for which Antioch default mixture file has untrusted \Delta H_f values
    void check_warn_for_species_with_untrusted_hf(const std::string& name) const;

        //! not allowed
        ASCIIParser();

        std::ifstream _doc;
        std::map<ParsingUnit,std::string> _unit_map;

        std::vector<unsigned int> _ignored;
        const unsigned int _n_columns_chemical_species;
        const unsigned int _n_columns_vib_data;
        const unsigned int _n_columns_el_data;
        const unsigned int _n_columns_transport_species;

    //! Needed to be able to warn about using not-trusted enthalpies of formation
    bool _is_antioch_default_mixture_file;
  };

} // end namespace Antioch


#endif
