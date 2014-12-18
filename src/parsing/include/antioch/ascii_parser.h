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
#include "antioch/parsing_enum.h"
#include "antioch/input_utils.h"
#include "antioch/string_utils.h"
#include "antioch/units.h"

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

  template <typename NumericType, typename CurveType>
  class NASAThermoMixture;

  // deprecated
  template <typename NumericType>
  class CEAThermodynamics;

  template <typename ThermoEvaluator, typename NumericType>
  class TransportMixture;

  template <typename NumericType>
  class ASCIIParser
  {
      public:
        ASCIIParser(const std::string& file, bool verbose = true);
        ~ASCIIParser();

        //! set the indexes of to-be-ignored columns
        void set_ignored_columns(const std::vector<unsigned int> & ignored);

        //! set the comment characters (overwrite the defaults: '#', '!')
        void set_comment_characters(const std::string & comments);

        //! read species list
        const std::vector<std::string>  species_list();

        //! read the mandatory data
        void read_chemical_species(ChemicalMixture<NumericType> & chem_mixture);

        //! read the vibrational data
        void read_vibrational_data(ChemicalMixture<NumericType>& chem_mixture);

        //! read the electronic data
        void read_electronic_data(ChemicalMixture<NumericType>& chem_mixture);

        //! read the thermodynamic data
        template <typename CurveType>
        void read_thermodynamic_data(NASAThermoMixture<NumericType, CurveType >& thermo);

        //! read the thermodynamic data, deprecated object
        void read_thermodynamic_data(CEAThermodynamics<NumericType>& thermo);

        //! read the transport data
        template <typename ThermoEvaluator>
        void read_transport_species(TransportMixture<ThermoEvaluator,NumericType> & transport_mixture);

        //! filename
        const std::string & filename() const;


     private:
        //! not allowed
        ASCIIParser();

        //! find the index of the wanted data
        void find_first(unsigned int & index,unsigned int n_data) const;

        //! several characters possible
        void skip_comments_lines();

        std::string   _file;
        std::ifstream _doc;
        bool          _verbose;
        std::map<ParsingUnit,std::string> _unit_map;
        std::vector<unsigned int>         _ignored;
        std::string                       _comments;

  };


  template <typename NumericType>
  inline
  ASCIIParser<NumericType>::ASCIIParser(const std::string & file, bool verbose):
        _file(file),
        _doc(file.c_str()),
        _verbose(verbose),
        _comments("#!")
  {
      if(!_doc.is_open())
      {
        std::cerr << "ERROR: unable to load file " << file << std::endl;
        antioch_error();
      }

      if(_verbose)std::cout << "Having opened file " << file << std::endl;

      skip_comments_lines();

      _unit_map[MOL_WEIGHT]    = "g/mol";
      _unit_map[MASS_ENTHALPY] = "J/kg";
  }

  template <typename NumericType>
  inline
  void ASCIIParser<NumericType>::skip_comments_lines()
  {
     for(unsigned int c = 0; c < _comments.size(); c++)
     {
        skip_comment_lines(_doc, _comments[c]);
     }
  }

  template <typename NumericType>
  inline
  void ASCIIParser<NumericType>::set_comment_characters(const std::string & comments)
  {
     _comments = comments;
  }

  template <typename NumericType>
  inline
  ASCIIParser<NumericType>::~ASCIIParser()
  {
     _doc.close();
  }

  template <typename NumericType>
  inline
  void ASCIIParser<NumericType>::set_ignored_columns(const std::vector<unsigned int> & ignored)
  {
      _ignored = ignored;
  }


  template <typename NumericType>
  inline
  const std::vector<std::string> ASCIIParser<NumericType>::species_list()
  {
      std::vector<std::string> species_list;
      std::string spec;

      while(_doc.good())
      {
          skip_comments_lines(); // if comments in the middle
          
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
          if(_verbose)std::cout << spec << std::endl;
          species_list.push_back(spec);
      }

      if(_verbose)std::cout << "Found " << species_list.size() << " species\n\n" << std::endl;
      return species_list;
  }

  template <typename NumericType>
  inline
  void ASCIIParser<NumericType>::find_first(unsigned int & index,unsigned int n_data) const
  {
        // slow algorithm so whatever order is fine
       bool find(true);
       while(find)
       {
          find = false;
          for(unsigned int ii = 0; ii < _ignored.size(); ii++)
          {
             if(index == _ignored[ii])
             {
                find = true;
                index++;
                break;
             }
          }
       }

       if(index > n_data)
       {
          std::cerr << "Error while reading " << _file << std::endl
                    << "Total number of columns provided is " << n_data
                    << " with " << _ignored.size() << " ignored column." << std::endl
                    << "The provided ignored index are:\n";
          for(unsigned int i = 0; i < _ignored.size(); i++)std::cerr << _ignored[i] << std::endl;
          std::cerr << "Indexes start at zero, maybe try decreasing them?" << std::endl;
          antioch_parsing_error("Error in ASCII parser");
       }
  }

  template <typename NumericType>
  inline
  void ASCIIParser<NumericType>::read_chemical_species(ChemicalMixture<NumericType> & chem_mixture)
  {
      
    std::string name;
    NumericType mol_wght, h_form, n_tr_dofs;
    int charge;
    NumericType mw_unit = Units<NumericType>(_unit_map.at(MOL_WEIGHT)).get_SI_factor();
    NumericType ef_unit = 1.L;// Units<NumericType>(_unit_map.at(MASS_ENTHALPY)).get_SI_factor(); // not integrated yet the kg bugfix

    const unsigned int n_data = 4 + _ignored.size(); // we read all those columns
    unsigned int imw(0);
    this->find_first(imw,n_data);
    unsigned int ihf(imw+1);
    this->find_first(ihf,n_data);
    unsigned int itrdofs(ihf+1);
    this->find_first(itrdofs,n_data);
    unsigned int ic(itrdofs+1);
    this->find_first(ic,n_data);

    std::vector<NumericType> read(n_data,0.);

    if(_verbose)std::cout << "Reading species characteristics in file " << _file << std::endl;
    while (_doc.good())
      {

         skip_comments_lines(); // if comment in the middle

        _doc >> name;      // Species Name
        for(unsigned int i = 0; i < n_data; i++)_doc >> read[i];
        mol_wght  = read[imw];     // molecular weight (kg/kmol)
        h_form    = read[ihf];     // heat of formation at Ok (J/kg)
        n_tr_dofs = read[itrdofs]; // number of translational/rotational DOFs
        charge    = int(read[ic]); // charge number

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
                if(_verbose)
                {
                    std::cout << "Adding " << name << " informations:\n\t"
                              << "molecular weight: "             << mol_wght  << " kg/mol\n\t"
                              << "formation enthalpy @0 K: "      << h_form    << " J/mol\n\t"
                              << "trans-rot degrees of freedom: " << n_tr_dofs << "\n\t"
                              << "charge: "                       << charge    << std::endl;
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

    const unsigned int n_data = 2 + _ignored.size(); // we read all those columns
    unsigned int itv(0);
    this->find_first(itv,n_data);
    unsigned int ide(itv+1);
    this->find_first(ide,n_data);

    std::vector<NumericType> read(n_data,0.);

    if(_verbose)std::cout << "Reading vibrational data in file " << _file << std::endl;
    while (_doc.good())
      {

        skip_comments_lines();

        _doc >> name;           // Species Name
        for(unsigned int i = 0; i < n_data; i++)_doc >> read[i];
        theta_v        = read[itv];                 // characteristic vibrational temperature (K)
        n_degeneracies = (unsigned int)(read[ide]); // degeneracy of the mode
      
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
            if(_verbose)
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

    const unsigned int n_data = 2 + _ignored.size(); // we read all those columns
    unsigned int ite(0);
    this->find_first(ite,n_data);
    unsigned int ide(ite+1);
    this->find_first(ide,n_data);

    std::vector<NumericType> read(n_data,0.);

    
    if(_verbose)std::cout << "Reading electronic data in file " << _file << std::endl;
    while (_doc.good())
      {
        _doc >> name;           // Species Name
        for(unsigned int i = 0; i < n_data; i++)_doc >> read[i];
        theta_e        = read[ite];                 // characteristic electronic temperature (K)
        n_degeneracies = (unsigned int)(read[ide]); // number of degeneracies for this level
        
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
            if(_verbose)
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
  void ASCIIParser<NumericType>::read_thermodynamic_data(NASAThermoMixture<NumericType, CurveType >& thermo)
  {
    std::string name;
    unsigned int n_int;
    std::vector<NumericType> coeffs;
    NumericType h_form, val;

    const ChemicalMixture<NumericType>& chem_mixture = thermo.chemical_mixture();

// \todo: only cea, should do NASA
    while (_doc.good())
      {
        skip_comments_lines();

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
                if(_verbose)std::cout << "Adding curve fit " << name << std::endl;
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
        skip_comments_lines();

        _doc >> name;   // Species Name
        _doc >> n_int;  // Number of T intervals: [200-1000], [1000-6000], ([6000-20000])
        _doc >> h_form; // Formation Enthalpy @ 298.15 K

        coeffs.clear();
        for (unsigned int interval=0; interval<n_int; interval++)
          {
            for (unsigned int n=0; n < 10; n++)
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
  template <typename ThermoEvaluator>
  inline
  void ASCIIParser<NumericType>::read_transport_species(TransportMixture<ThermoEvaluator,NumericType> & transport_mixture)
  {
      
    std::string name;
    NumericType LJ_eps_kB;
    NumericType LJ_sigma;
    NumericType dipole_moment;
    NumericType pol;
    NumericType Zrot;

    const unsigned int n_data = 5 + _ignored.size(); // we read all those columns
    unsigned int iLJeps(0);
    this->find_first(iLJeps,n_data);
    unsigned int iLJsig(iLJeps+1);
    this->find_first(iLJsig,n_data);
    unsigned int idip(iLJsig+1);
    this->find_first(idip,n_data);
    unsigned int ipol(idip+1);
    this->find_first(ipol,n_data);
    unsigned int irot(ipol+1);
    this->find_first(irot,n_data);

    std::vector<NumericType> read(n_data,0.);

    while (_doc.good())
    {
        skip_comments_lines(); // if comment in the middle
        _doc >> name;
        for(unsigned int i = 0; i < n_data; i++)_doc >> read[i];
        LJ_eps_kB     = read[iLJeps];
        LJ_sigma      = read[iLJsig];
        dipole_moment = read[idip];
        pol           = read[ipol];
        Zrot          = read[irot];

        if(transport_mixture.chemical_mixture().species_name_map().count(name))
        {
           unsigned int place = transport_mixture.chemical_mixture().species_name_map().at(name);

           NumericType mass = transport_mixture.chemical_mixture().M(place);
           // adding species in mixture
           transport_mixture.add_species(place,LJ_eps_kB,LJ_sigma,dipole_moment,pol,Zrot,mass);
        }else
        {
        }
    }
  }

  template <typename NumericType>
  inline
  const std::string & ASCIIParser<NumericType>::filename() const
  {
     return _file;
  }
  

} // end namespace Antioch


#endif
