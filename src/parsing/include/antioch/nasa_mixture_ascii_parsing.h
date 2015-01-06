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
//--------------------------------------------------------------------------

#ifndef ANTIOCH_NASA_MIXTURE_ASCII_PARSING_H
#define ANTIOCH_NASA_MIXTURE_ASCII_PARSING_H

// C++
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

// Antioch
#include "antioch/input_utils.h"
#include "antioch/chemical_mixture.h"

namespace Antioch
{
  // Forward declarations
  template <class NumericType, class NASAFit>
  class NASAThermoMixture;

  template <class NumericType>
  class NASA7CurveFit;

  // New declarations

  template<class NumericType>
  void read_nasa_mixture_data_ascii( NASAThermoMixture<NumericType,NASA7CurveFit<NumericType> >& thermo, const std::string &filename );

 
  /* ------------------------- Inline Functions -------------------------*/
  template<class NumericType>
  inline
  void read_nasa_mixture_data_ascii( NASAThermoMixture<NumericType, NASA7CurveFit<NumericType> >& thermo, const std::string &filename )
  {
    
// ChemKin style, except no headers, direct to
// the species, chemkin classic
// only two intervals
// TODO: the parser will allow custom
// intervals definition
    std::ifstream in(filename.c_str());
    if(!in.is_open())
    {
      std::cerr << "ERROR: unable to load file " << filename << std::endl;
      antioch_error();
    }

    skip_comment_lines(in, '!');

    std::string name;
    std::vector<NumericType> coeffs;
    std::vector<NumericType> temps(3,0.);

    const ChemicalMixture<NumericType>& chem_mixture = thermo.chemical_mixture();

    while (in.good())
      {
        std::string line;
        std::stringstream tmp;

        if(!getline(in,line))break;
        if(line[0] == '!')continue;
        tmp << line.substr(0,18); //name

        // this is ChemKin doc, 1st: temps E10.0
        tmp << line.substr(45,10);         // low
        tmp << " " << line.substr(55,10);  // high
        tmp << " " << line.substr(65,10);  // inter

         // get rid of last character
        for(unsigned int n = 0; n < 3; n++)
        {
          if(!getline(in,line))// we have a problem
          {
             std::cerr << "NASA input file error, check " << filename << std::endl;
             in.close();
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
	if (in.good())
	  {
	    // Check if this is a species we want.
	    if( chem_mixture.species_name_map().find(name) !=
		chem_mixture.species_name_map().end() )
	      {
		thermo.add_curve_fit(name, coeffs, temps);
	      }
	  }
      } // end while
 
   
    // Make sure we actually populated everything
    if( !thermo.check() )
      {
	std::cerr << "Error: NASA table not fully populated" << std::endl;
	antioch_error();
      }

    in.close();

    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_NASA_MIXTURE_ASCII_PARSING_H
