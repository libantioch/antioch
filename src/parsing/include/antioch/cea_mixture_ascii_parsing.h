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

#ifndef ANTIOCH_NASA9_MIXTURE_ASCII_PARSING_H
#define ANTIOCH_NASA9_MIXTURE_ASCII_PARSING_H

// C++
#include <iostream>
#include <vector>

// Antioch
#include "antioch/input_utils.h"
#include "antioch/chemical_mixture.h"

namespace Antioch
{
  // Forward declarations
  template <class NumericType, class NASAFit>
  class NASAThermoMixture;

  template <class NumericType>
  class NASA9CurveFit;

  // New declarations

  template<class NumericType>
  void read_cea_mixture_data_ascii( NASAThermoMixture<NumericType,NASA9CurveFit<NumericType> >& thermo, const std::string &filename );

  template<class NumericType>
  void read_cea_mixture_data_ascii_default( NASAThermoMixture<NumericType,NASA9CurveFit<NumericType> >& thermo );

 
  /* ------------------------- Inline Functions -------------------------*/
  template<class NumericType>
  inline
  void read_cea_mixture_data_ascii( NASAThermoMixture<NumericType, NASA9CurveFit<NumericType> >& thermo, const std::string &filename )
  {
    
    std::ifstream in(filename.c_str());
    if(!in.is_open())
    {
      std::cerr << "ERROR: unable to load file " << filename << std::endl;
      antioch_error();
    }

    skip_comment_lines(in, '#');

    std::string name;
    unsigned int n_int;
    std::vector<NumericType> coeffs;
    NumericType h_form, val;

    const ChemicalMixture<NumericType>& chem_mixture = thermo.chemical_mixture();

    while (in.good())
      {
	in >> name;   // Species Name
	in >> n_int;  // Number of T intervals: [200-1000], [1000-6000], ([6000-20000])
	in >> h_form; // Formation Enthalpy @ 298.15 K

	coeffs.clear();
	for (unsigned int interval=0; interval<n_int; interval++)
	  {
	    for (unsigned int n=0; n<10; n++)
	      {
		in >> val, coeffs.push_back(val);
	      }
	  }

	// If we are still good, we have a valid set of thermodynamic
	// data for this species. Otherwise, we read past end-of-file 
	// in the section above
	if (in.good())
	  {
	    // Check if this is a species we want.
	    if( chem_mixture.species_name_map().find(name) !=
		chem_mixture.species_name_map().end() )
	      {
		thermo.add_curve_fit(name, coeffs);
	      }
	  }
      } // end while
    
    // Make sure we actually populated everything
    if( !thermo.check() )
      {
	std::cerr << "Error: NASA9 table not fully populated" << std::endl;
	antioch_error();
      }

    in.close();

    return;
  }

  template<class NumericType>
  void read_cea_mixture_data_ascii_default( NASAThermoMixture<NumericType, NASA9CurveFit<NumericType> >& thermo )
  {
    antioch_deprecated();
    read_cea_mixture_data_ascii(thermo, DefaultFilename::thermo_data());
  }

} // end namespace Antioch

#endif // ANTIOCH_NASA9_MIXTURE_ASCII_PARSING_H
