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

#ifndef ANTIOCH_ASCII_PARSER_H
#define ANTIOCH_ASCII_PARSER_H

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/parser_base.h"
#include "antioch/parsing_enum.h"
#include "antioch/input_utils.h"
#include "antioch/string_utils.h"
#include "antioch/units.h"
#include "antioch/cea_thermo.h" // because not templated, therefore should be entirely known
#include "antioch/cea_curve_fit.h" // because not templated, therefore should be entirely known
#include "antioch/nasa_curve_fit.h" // because not templated, therefore should be entirely known
#include "antioch/transport_mixture.h" // because not templated, therefore should be entirely known

// C++
#include <iostream>
#include <fstream>
#include <algorithm> // std::search_n
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

  template <class ThermoEval,class NumericType>
  class TransportMixture;

  template <typename NumericType, typename Macro,typename Micro>
  class ThermoHandler;

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
  class CEAThermodynamics;

  template <typename NumericType>
  class CEAEvaluator;

// micro
  template <typename NumericType>
  class StatMechThermodynamics;

  template <typename Macro, typename NumericType>
  class IdealGasKineticsTheory;

  template <typename NumericType>
  class ASCIIParser: public ParserBase<NumericType>
  {
      public:
        ASCIIParser(const std::string& file, bool verbose = true);
        ~ASCIIParser();

        void change_file(const std::string & filename);

        bool initialize() {return false;}

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

// transport, the thermo is explicit...

        //! reads the transport data, not valid in xml && chemkin
        //  NASA7 + StatMech
        void read_transport_data(TransportMixture< ThermoHandler < NumericType,
                                                                   NASAEvaluator<NumericType,NASA7CurveFit<NumericType> >,
                                                                   StatMechThermodynamics<NumericType> 
                                                                 > ,
                                                   NumericType > & transport_mixture) {this->read_transport_data_root(transport_mixture);}

        //! reads the transport data, not valid in xml && chemkin
        //  NASA9 + StatMech
        void read_transport_data(TransportMixture< ThermoHandler < NumericType,
                                                                   NASAEvaluator<NumericType,NASA9CurveFit<NumericType> >,
                                                                   StatMechThermodynamics<NumericType>
                                                                 >,
                                                   NumericType > & transport_mixture) {this->read_transport_data_root(transport_mixture);}

        //! reads the transport data, not valid in xml && chemkin
        //  CEA + StatMech for backward compat
        void read_transport_data(TransportMixture< ThermoHandler < NumericType,
                                                                   CEAEvaluator<NumericType>,
                                                                   StatMechThermodynamics<NumericType> 
                                                                 >,
                                                   NumericType > & transport_mixture) {this->read_transport_data_root(transport_mixture);}

        //! reads the transport data, not valid in xml && chemkin
        //  NASA7 + Ideal Gas
        void read_transport_data(TransportMixture< ThermoHandler < NumericType,
                                                                   NASAEvaluator<NumericType,NASA7CurveFit<NumericType> >,
                                                                   IdealGasKineticsTheory<NASAEvaluator<NumericType,NASA7CurveFit<NumericType> >, NumericType> 
                                                                 >,
                                                   NumericType > & transport_mixture) {this->read_transport_data_root(transport_mixture);}

        //! reads the transport data, not valid in xml && chemkin
        //  NASA9 + Ideal Gas
        void read_transport_data(TransportMixture< ThermoHandler < NumericType,
                                                                   NASAEvaluator<NumericType,NASA9CurveFit<NumericType> >,
                                                                   IdealGasKineticsTheory<NASAEvaluator<NumericType,NASA9CurveFit<NumericType> >, NumericType> 
                                                                 >,
                                                   NumericType > & transport_mixture) {this->read_transport_data_root(transport_mixture);}

        //! reads the transport data, not valid in xml && chemkin
        //  CEA + Ideal Gas for backward compat
        void read_transport_data(TransportMixture< ThermoHandler < NumericType,
                                                                   CEAEvaluator<NumericType>,
                                                                   IdealGasKineticsTheory<CEAEvaluator<NumericType>,NumericType>
                                                                 >,
                                                   NumericType > & transport_mixture) {this->read_transport_data_root(transport_mixture);}


///////////////
// thermo
///////////////

//global overload
        //! reads the thermo, NASA generalist, no templates for virtual
        void read_thermodynamic_data(NASAThermoMixture<NumericType, NASA7CurveFit<NumericType> >& thermo)
                {this->read_thermodynamic_data_root(thermo);}

        //! reads the thermo, NASA generalist, no templates for virtual
        void read_thermodynamic_data(NASAThermoMixture<NumericType, NASA9CurveFit<NumericType> >& thermo)
                {this->read_thermodynamic_data_root(thermo);}

        //! reads the thermo, NASA generalist, no templates for virtual
        void read_thermodynamic_data(NASAThermoMixture<NumericType, CEACurveFit<NumericType> >& thermo)  
                {this->read_thermodynamic_data_root(thermo);}

        //! read the thermodynamic data, deprecated object
        void read_thermodynamic_data(CEAThermodynamics<NumericType>& thermo);



//////////////
//kinetics
//////////////


     private:

// templated thermo version
        //! read the thermodynamic data
        template <typename CurveType>
        void read_thermodynamic_data_root(NASAThermoMixture<NumericType, CurveType >& thermo);

// templated transport version
        //! read the thermodynamic data
        template <typename Mixture>
        void read_transport_data_root(Mixture & transport);

        //! not allowed
        ASCIIParser();

        std::ifstream _doc;
        std::map<ParsingUnit,std::string> _unit_map;

  };


  template <typename NumericType>
  inline
  ASCIIParser<NumericType>::ASCIIParser(const std::string & file, bool verbose):
        ParserBase<NumericType>("ascii",file,verbose),
        _doc(file.c_str())
  {
      if(!_doc.is_open())
      {
        std::cerr << "ERROR: unable to load file " << file << std::endl;
        antioch_error();
      }

      if(this->verbose())std::cout << "Having opened file " << file << std::endl;

      skip_comment_lines(_doc, '#');

      _unit_map[MOL_WEIGHT]    = "g/mol";
      _unit_map[MASS_ENTHALPY] = "J/kg";
  }

  template <typename NumericType>
  inline
  ASCIIParser<NumericType>::~ASCIIParser()
  {
     _doc.close();
  }

  template <typename NumericType>
  inline
  void ASCIIParser<NumericType>::change_file(const std::string & filename)
  {
    _doc.close();
    _doc.open(filename.c_str());
    ParserBase<NumericType>::_file = filename;
    if(!_doc.is_open())
      {
        std::cerr << "ERROR: unable to load file " << filename << std::endl;
        antioch_error();
      }

      if(this->verbose())std::cout << "Having opened file " << filename << std::endl;

      skip_comment_lines(_doc, '#');
  }


  template <typename NumericType>
  inline
  const std::vector<std::string> ASCIIParser<NumericType>::species_list()
  {
      std::vector<std::string> species_list;
      std::string spec;

      while(_doc.good())
      {
          skip_comment_lines(_doc, '#'); // if comments in the middle
          
          _doc >> spec;

          if(!_doc.good())break; // read successful?

          // checking doublon
          bool doublon(false);
          for(unsigned int s = 0; s < species_list.size(); s++)
          {
              if(spec == species_list[s])
              {
                doublon = true;
                break;
              }
              if(doublon)
              {
                 std::cerr << "Multiple declaration of " << spec
                           << ", skipping doublon" << std::endl;
                 continue;
              }
              // adding
          }
          if(this->verbose())std::cout << spec << std::endl;
          species_list.push_back(spec);
      }

      if(this->verbose())std::cout << "Found " << species_list.size() << " species\n\n" << std::endl;
      return species_list;
  }

  template <typename NumericType>
  inline
  void ASCIIParser<NumericType>::read_chemical_species(ChemicalMixture<NumericType> & chem_mixture)
  {
      
    std::string name;
    NumericType mol_wght, h_form, n_tr_dofs;
    int charge;
    NumericType mw_unit = Units<NumericType>(_unit_map.at(MOL_WEIGHT)).get_SI_factor();
    NumericType ef_unit = Units<NumericType>(_unit_map.at(MASS_ENTHALPY)).get_SI_factor(); // not integrated yet the kg bugfix

    if(this->verbose())std::cout << "Reading species characteristics in file " << this->file() << std::endl;
    while (_doc.good())
      {

         skip_comment_lines(_doc, '#'); // if comment in the middle

        _doc >> name;      // Species Name
        _doc >> mol_wght;  // molecular weight (kg/kmol)
        _doc >> h_form;    // heat of formation at Ok (J/kg)
        _doc >> n_tr_dofs; // number of translational/rotational DOFs
        _doc >> charge;    // charge number

        mol_wght *= mw_unit; // to SI (kg/mol)
        h_form *= ef_unit; // to SI (J/kg)
        
        // If we are still good, we have a valid set of thermodynamic
        // data for this species. Otherwise, we read past end-of-file 
        // in the section above
        if (_doc.good())
          {
            // If we do not have this species, just go on
            if (!chem_mixture.species_name_map().count(name))continue;

            Species species = chem_mixture.species_name_map().at(name);

            // using default comparison:
            std::vector<Species>::const_iterator it = std::search_n( chem_mixture.species_list().begin(), 
                                                                     chem_mixture.species_list().end(),
                                                                     1, species);
            if( it != chem_mixture.species_list().end() )
              {
                unsigned int index = static_cast<unsigned int>(it - chem_mixture.species_list().begin());
                chem_mixture.add_species( index, name, mol_wght, h_form, n_tr_dofs, charge );
                if(this->verbose())
                {
                    std::cout << "Adding " << name << " informations:\n\t"
                              << "molecular weight: "             << mol_wght << " kg/mol\n\t"
                              << "formation enthalpy @0 K: "      << h_form << " J/mol\n\t"
                              << "trans-rot degrees of freedom: " << n_tr_dofs << "\n\t"
                              << "charge: "                       << charge << std::endl;
                }
              }

          }
      }

      return;
  }

  template <typename NumericType>
  inline
  void ASCIIParser<NumericType>::read_vibrational_data(ChemicalMixture<NumericType> & chem_mixture)
  {
    std::string name;
    NumericType theta_v;
    unsigned int n_degeneracies;

    if(this->verbose())std::cout << "Reading vibrational data in file " << this->file() << std::endl;
    while (_doc.good())
      {

        skip_comment_lines(_doc, '#'); // if comment in the middle

        _doc >> name;           // Species Name
        _doc >> theta_v;        // characteristic vibrational temperature (K)
        _doc >> n_degeneracies; // degeneracy of the mode
      
        // If we are still good, we have a valid set of thermodynamic
        // data for this species. Otherwise, we read past end-of-file 
        // in the section above
        if (_doc.good())
          {
            // If we do not have this species, just keep going
            if (!chem_mixture.species_name_map().count(name))continue;
        
            // ... otherwise we add the data
            const unsigned int s = 
              chem_mixture.species_name_map().find(name)->second;

            antioch_assert_equal_to((chem_mixture.chemical_species()[s])->species(), name);
        
            chem_mixture.add_species_vibrational_data(s, theta_v, n_degeneracies);
            if(this->verbose())
            {
                std::cout << "Adding vibrational data of species " << name << "\n\t"
                          << "vibrational temperature: " << theta_v << " K\n\t"
                          << "degeneracy: "              << n_degeneracies << std::endl;
            }
          }
      }

      return;
  }

  template <typename NumericType>
  inline
  void ASCIIParser<NumericType>::read_electronic_data(ChemicalMixture<NumericType> & chem_mixture)
  {
    std::string name;
    NumericType theta_e;
    unsigned int n_degeneracies;
    
    if(this->verbose())std::cout << "Reading electronic data in file " << this->file() << std::endl;
    while (_doc.good())
      {
        _doc >> name;           // Species Name
        _doc >> theta_e;        // characteristic electronic temperature (K)
        _doc >> n_degeneracies; // number of degeneracies for this level
        
        // If we are still good, we have a valid set of thermodynamic
        // data for this species. Otherwise, we read past end-of-file 
        // in the section above
        if (_doc.good())
          {
            // If we do not have this species, just go on
            if (!chem_mixture.species_name_map().count(name))continue;
            
            // ... otherwise we add the data
            const unsigned int s = 
              chem_mixture.species_name_map().find(name)->second;

            antioch_assert_equal_to((chem_mixture.chemical_species()[s])->species(), name);
            
            chem_mixture.add_species_electronic_data(s, theta_e, n_degeneracies);
            if(this->verbose())
            {
                std::cout << "Adding electronic data of species " << name << "\n\t"
                          << "electronic temperature: " << theta_e << " K\n\t"
                          << "degeneracy: "             << n_degeneracies << std::endl;
            }
          }
      }
      return;
  }

  template <typename NumericType>
  template <typename CurveType>
  inline
  void ASCIIParser<NumericType>::read_thermodynamic_data_root(NASAThermoMixture<NumericType, CurveType >& thermo)
  {
    std::string name;
    unsigned int n_int;
    std::vector<NumericType> coeffs;
    NumericType h_form, val;

    const ChemicalMixture<NumericType>& chem_mixture = thermo.chemical_mixture();

// \todo: only cea, should do NASA
    while (_doc.good())
      {
        skip_comment_lines(_doc, '#');

        _doc >> name;   // Species Name
        _doc >> n_int;  // Number of T intervals: [200-1000], [1000-6000], ([6000-20000])
        _doc >> h_form; // Formation Enthalpy @ 298.15 K

        coeffs.clear();
        for (unsigned int interval=0; interval<n_int; interval++)
          {
            for (unsigned int n=0; n<10; n++)
              {
                _doc >> val, coeffs.push_back(val);
              }
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
                if(this->verbose())std::cout << "Adding curve fit " << name << std::endl;
                thermo.add_curve_fit(name, coeffs);
              }
          }
      } // end while
  }


  template <typename NumericType>
  inline
  void ASCIIParser<NumericType>::read_thermodynamic_data(CEAThermodynamics<NumericType>& thermo)
  {
    std::string name;
    unsigned int n_int;
    std::vector<NumericType> coeffs;
    NumericType h_form, val;

    const ChemicalMixture<NumericType>& chem_mixture = thermo.chemical_mixture();

// \todo: only cea, should do NASA
    while (_doc.good())
      {
        skip_comment_lines(_doc, '#');

        _doc >> name;   // Species Name
        _doc >> n_int;  // Number of T intervals: [200-1000], [1000-6000], ([6000-20000])
        _doc >> h_form; // Formation Enthalpy @ 298.15 K

        coeffs.clear();
        for (unsigned int interval=0; interval<n_int; interval++)
          {
            for (unsigned int n=0; n<10; n++)
              {
                _doc >> val, coeffs.push_back(val);
              }
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
                thermo.add_curve_fit(name, coeffs);
              }
          }
      } // end while
  }

  template <typename NumericType>
  template <typename Mixture>
  inline
  void ASCIIParser<NumericType>::read_transport_data_root(Mixture & transport)
  {
        
    std::string name;
    NumericType LJ_eps_kB;
    NumericType LJ_sigma;
    NumericType dipole_moment;
    NumericType pol;
    NumericType Zrot;
    NumericType SI_coeff = Antioch::Units<NumericType>("g/mol").get_SI_factor();

    while (_doc.good())
      {
          skip_comment_lines(_doc, '#');

          _doc >> name >> LJ_eps_kB >> LJ_sigma >> dipole_moment>> pol >> Zrot;
          if(transport.chemical_mixture().species_name_map().count(name))
          {
              unsigned int place = transport.chemical_mixture().species_name_map().at(name);
              NumericType mass = transport.chemical_mixture().M(place) * SI_coeff;
// adding species in mixture
              transport.add_species(place,LJ_eps_kB,LJ_sigma,dipole_moment,pol,Zrot,mass);
          }
      }
  }

} // end namespace Antioch


#endif
