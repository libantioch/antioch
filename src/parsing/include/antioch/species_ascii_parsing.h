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

#ifndef ANTIOCH_SPECIES_ASCII_PARSING_H
#define ANTIOCH_SPECIES_ASCII_PARSING_H

// Antioch
#include "antioch/input_utils.h"
#include "antioch/string_utils.h"

// C++
#include <algorithm> // std::search_n
#include <fstream>

namespace Antioch
{
  typedef unsigned int Species;

  // Forward declarations
  template <class NumericType>
  class ChemicalMixture;

  template <typename NumericType>
  void read_chemical_species_composition(const std::string & filename,
                                         bool verbose,
                                         ChemicalMixture<NumericType> & mixture);

  template<class NumericType>
  void read_species_data_ascii( const std::string& filename,
                                bool verbose,
                                ChemicalMixture<NumericType>& chem_mixture);
				

  template<class NumericType>
  void read_species_vibrational_data_ascii (const std::string &filename,
                                            bool verbose,
                                            ChemicalMixture<NumericType>& chem_mixture);

  template<class NumericType>
  void read_species_electronic_data_ascii (const std::string &filename,
                                           bool verbose,
                                           ChemicalMixture<NumericType>& chem_mixture);

//----------------------------------------------------------------

   template <typename NumericType>
   inline
   void read_chemical_species_composition(const std::string & filename,
                                          bool verbose,
                                          ChemicalMixture<NumericType> & mixture)
   {
      std::ifstream data(filename.c_str());
      if(!data.is_open())
      {
        std::cerr << "ERROR: unable to load file " << filename << std::endl;
        antioch_error();
      }

      skip_comment_lines(data, '#');

      std::vector<std::string> species_list;
      std::string line;
      while(!data.eof())
      {
          if(!getline(data,line))break;
          if(line[0] == '#' || line.empty())continue;
          std::vector<std::string> tmp_spec;
// in case several species on same line
          int nspec = SplitString(line," ",tmp_spec,false);
// in case only one
          if(nspec == 0)tmp_spec.push_back(line);
          for(unsigned int i = 0; i < tmp_spec.size(); i++)
          {
              // checking if multiple declaration
              bool doublon(false);
              for(unsigned int s = 0; s < species_list.size(); s++)
              {
                 if(tmp_spec[i] == species_list[s])
                 {
                    doublon = true;
                    break;
                 }
              }
              if(doublon)
              {
                 std::cerr << "Multiple declaration of " << tmp_spec[i]
                           << ", skipping doublon" << std::endl;
                 continue;
              }
              // adding
              if(verbose)std::cout << tmp_spec[i] << std::endl;
              species_list.push_back(tmp_spec[i]);
          }
      }
      data.close();

      if(verbose)std::cout << "Added " << species_list.size() << " species" << std::endl;

      mixture.initialize_species(species_list);
   }
  
  template<class NumericType>
  inline
  void read_species_data_ascii(const std::string& filename, bool verbose, 
                               ChemicalMixture<NumericType>& chem_mixture)
  {
    std::ifstream in(filename.c_str());
    if(!in.is_open())
    {
       std::cerr << "ERROR: unable to load file " << filename << std::endl;
       antioch_error();
    }
    
    skip_comment_lines(in, '#');

    std::string name;
    NumericType mol_wght, h_form, n_tr_dofs;
    int charge;

    while (in.good())
      {
	in >> name;      // Species Name
	in >> mol_wght;  // molecular weight (kg/kmol)
	in >> h_form;    // heat of formation at Ok (J/kg)
	in >> n_tr_dofs; // number of translational/rotational DOFs
	in >> charge;    // charge number

        mol_wght *= 1e-3L; // to SI (kg/mol)
	
	// If we are still good, we have a valid set of thermodynamic
	// data for this species. Otherwise, we read past end-of-file 
	// in the section above
	if (in.good())
	  {
	    // If we do not have this species, just go on
	    if (!chem_mixture.species_name_map().count(name))continue;


	    /* Only add a ChemicalSpecies if the user supplied it to
	       the ChemicalMixture. */
	    Species species = chem_mixture.species_name_map().find(name)->second;

	    // using default comparison:
	    std::vector<Species>::const_iterator it = std::search_n( chem_mixture.species_list().begin(), 
								     chem_mixture.species_list().end(),
								     1, species);
	    if( it != chem_mixture.species_list().end() )
	      {
		unsigned int index = static_cast<unsigned int>(it - chem_mixture.species_list().begin());
		chem_mixture.add_species( index, name, mol_wght, h_form, n_tr_dofs, charge );
                if(verbose)
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
    in.close();

    // sanity check, we require these informations
    bool fail(false);
    for(unsigned int s = 0; s < chem_mixture.chemical_species().size(); s++)
    {
        if(!chem_mixture.chemical_species()[s])
        {
            fail = true;
            break;
        }
    }
    if(fail)
    {
      std::cerr << "Molecule(s) is(are) missing.  Please update the information."
                << "  Currently using file " << filename << ".\n"
                << "Missing molecule(s) is(are):" << std::endl;
      for(unsigned int i = 0; i < chem_mixture.species_list().size(); i++)
      {
        if(!chem_mixture.chemical_species()[i])
        {
           std::cerr << chem_mixture.species_inverse_name_map().at(i) << std::endl;
        }
      }
      antioch_error();
    }

    return;
  }

  template<class NumericType>
  inline
  void read_species_vibrational_data_ascii (const std::string & filename, 
                                            bool verbose, 
                                            ChemicalMixture<NumericType>& chem_mixture)
  {
    std::ifstream in(filename.c_str());
    if(!in.is_open())
    {
       std::cerr << "ERROR: unable to load file " << filename << std::endl;
       antioch_error();
    }
    // skip the header
    skip_comment_lines(in, '#');

    std::string name;
    NumericType theta_v;
    unsigned int n_degeneracies;

    while (in.good())
      {
        in >> name;           // Species Name
        in >> theta_v;        // characteristic vibrational temperature (K)
        in >> n_degeneracies; // degeneracy of the mode
      
        // If we are still good, we have a valid set of thermodynamic
        // data for this species. Otherwise, we read past end-of-file 
        // in the section above
        if (in.good())
          {
            // If we do not have this species, just keep going
            if (!chem_mixture.species_name_map().count(name))continue;
	
            // ... otherwise we add the data
            const unsigned int s = 
              chem_mixture.species_name_map().find(name)->second;

            antioch_assert_equal_to((chem_mixture.chemical_species()[s])->species(), name);
        
            chem_mixture.add_species_vibrational_data(s, theta_v, n_degeneracies);
            if(verbose)
            {
                std::cout << "Adding vibrational data of species " << name << "\n\t"
                          << "vibrational temperature: " << theta_v << " K\n\t"
                          << "degeneracy: "              << n_degeneracies << std::endl;
            }
          }
      }
      in.close();

    // sanity check, we check these informations
    std::vector<std::string> missing;
    for(unsigned int s = 0; s < chem_mixture.chemical_species().size(); s++)
    {
        if(chem_mixture.chemical_species()[s]->theta_v().empty())missing.push_back(chem_mixture.chemical_species()[s]->species());
    }
    if(!missing.empty())
    {
       std::cerr << "WARNING:\nVibrational levels are missing.  Please update the information."
                 << "  Currently using file " << filename << ".\n"
                 << "Missing molecule(s) is(are):" << std::endl;
       for(unsigned int m = 0; m < missing.size(); m++)std::cerr << missing[m] << std::endl;
    }
    return;
  }

  
  template<class NumericType>
  inline
  void read_species_electronic_data_ascii (const std::string & filename, bool verbose,
                                           ChemicalMixture<NumericType>& chem_mixture)
                                           
  {
    std::ifstream in(filename.c_str());
    if(!in.is_open())
    {
       std::cerr << "ERROR: unable to load file " << filename << std::endl;
       antioch_error();
    }
    // skip the header
    skip_comment_lines(in, '#');
    
    std::string name;
    NumericType theta_e;
    unsigned int n_degeneracies;
    
    while (in.good())
      {
        in >> name;           // Species Name
        in >> theta_e;        // characteristic electronic temperature (K)
        in >> n_degeneracies; // number of degeneracies for this level
        
        // If we are still good, we have a valid set of thermodynamic
        // data for this species. Otherwise, we read past end-of-file 
        // in the section above
        if (in.good())
          {
            // If we do not have this species, just go on
            if (!chem_mixture.species_name_map().count(name))continue;
            
            // ... otherwise we add the data
            const unsigned int s = 
              chem_mixture.species_name_map().find(name)->second;

            antioch_assert_equal_to((chem_mixture.chemical_species()[s])->species(), name);
            
            chem_mixture.add_species_electronic_data(s, theta_e, n_degeneracies);
            if(verbose)
            {
                std::cout << "Adding electronic data of species " << name << "\n\t"
                          << "electronic temperature: " << theta_e << " K\n\t"
                          << "degeneracy: "             << n_degeneracies << std::endl;
            }
          }
        
      }
      in.close();
    // sanity check, we check these informations
    std::vector<std::string> missing;
    for(unsigned int s = 0; s < chem_mixture.chemical_species().size(); s++)
    {
        if(chem_mixture.chemical_species()[s]->theta_e().empty())missing.push_back(chem_mixture.chemical_species()[s]->species());
    }
    if(!missing.empty())
    {
       std::cerr << "WARNING:\nElectronic levels are missing.  Please update the information."
                 << "  Currently using file " << filename << ".\n"
                 << "Missing molecule(s) is(are):" << std::endl;
       for(unsigned int m = 0; m < missing.size(); m++)std::cerr << missing[m] << std::endl;
    }
    return;
  }
  
} // end namespace Antioch

#endif // ANTIOCH_SPECIES_ASCII_PARSING_H
