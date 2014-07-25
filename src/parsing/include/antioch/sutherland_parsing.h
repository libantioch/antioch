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

#ifndef ANTIOCH_SUTHERLAND_PARSING_H
#define ANTIOCH_SUTHERLAND_PARSING_H

// Antioch
#include "antioch/sutherland_viscosity.h"
#include "antioch/mixture_viscosity.h"

// C++
#include <iostream>

namespace Antioch
{
  template<class NumericType>
  void read_sutherland_data_ascii( MixtureViscosity<SutherlandViscosity<NumericType>,NumericType >& mu,
				   std::istream &in);

  template<class NumericType>
  void read_sutherland_data_ascii_default( MixtureViscosity<SutherlandViscosity<NumericType>,NumericType>& mu );

 
  /* ------------------------- Inline Functions -------------------------*/
  template<class NumericType>
  inline
  void read_sutherland_data_ascii( MixtureViscosity<SutherlandViscosity<NumericType>,NumericType >& mu,
				   std::istream &in)
  {
    // skip the header
    skip_comment_lines(in, '#');

    std::string name;
    NumericType a, b;

     while (in.good())
      {
        in >> name; // Species Name
        in >> a;    // 
        in >> b;    //
        
        // If we are still good, we have a valid set of transport
        // data for this species. Otherwise, we read past end-of-file 
        // in the section above
        if (in.good())
          {
	    const ChemicalMixture<NumericType>& chem_mixture = mu.chemical_mixture();
	    
	    // Check if this is a species we want.
	    if( chem_mixture.active_species_name_map().find(name) !=
		chem_mixture.active_species_name_map().end() )
	      {
		// Pack up coefficients
		std::vector<NumericType> coeffs(2);
		coeffs[0] = a;
		coeffs[1] = b;
		mu.add(name, coeffs);
	      }
	  }
      }

    return;
  }

  template<class NumericType>
  inline
  void read_sutherland_data_ascii_default( MixtureViscosity<SutherlandViscosity<NumericType>,NumericType>& mu )
  {
    static const std::string
      default_sutherland_transport_data
      ("#-----------------------------------------------------------------------------\n"
       "# Coefficients for Sutherland viscosity model\n"
       "# \n"
       "# Form of the fit:\n"
       "# \n"
       "# mu = A*T^1.5/(T + B)\n"
       "# \n"
       "# where T is in Kelvin and the viscosity is then given in Pa-s.\n"
       "# \n"
       "# Sources:\n"
       "# \n"
       "# Air   -- \n"
       "# N2    -- \n"
       "# CPAir & CPN2 -- compatibility with benkirk's dissertation cases\n"
       "# \n"
       "# Sutherland coefficients are in general not as accurate as Blottner fits,\n"
       "# and are provided mainly for completeness.\n"
       "\n"
       "Air     1.458000e-06  1.103000e+02\n"
       "CPAir   1.458000e-06  1.104000e+02\n"
       "N2      1.399306e-06  1.066667e+02\n"
       "CPN2    1.399306e-06  1.066667e+02\n");

    std::istringstream buf(default_sutherland_transport_data);

    read_sutherland_data_ascii( mu, buf );

    return;
  }


} // end namespace Antioch

#endif // ANTIOCH_SUTHERLAND_PARSING_H
