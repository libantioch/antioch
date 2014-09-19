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
#include "antioch/chemical_mixture.h"
#include "antioch/physical_set.h"

// C++
#include <iostream>

namespace Antioch
{
  template<class NumericType>
  void read_sutherland_data_ascii( PhysicalSet<SutherlandViscosity<NumericType>, ChemicalMixture<NumericType> >& mu,
				   const std::string &filename);
 
  template<class NumericType>
  void read_sutherland_data_ascii_default( PhysicalSet<SutherlandViscosity<NumericType>, ChemicalMixture<NumericType> >& mu );

  /* ------------------------- Inline Functions -------------------------*/
  template<class NumericType>
  inline
  void read_sutherland_data_ascii( PhysicalSet<SutherlandViscosity<NumericType>, ChemicalMixture<NumericType> >& mu,
				   const std::string &filename)
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
	    const ChemicalMixture<NumericType>& chem_mixture = mu.mixture();
	    
	    // Check if this is a species we want.
	    if( chem_mixture.species_name_map().find(name) !=
		chem_mixture.species_name_map().end() )
	      {
		// Pack up coefficients
		std::vector<NumericType> coeffs(2);
		coeffs[0] = a;
		coeffs[1] = b;
		mu.add_model(name, coeffs);
	      }
	  }
      }
      in.close();
    return;
  }

  template<class NumericType>
  void read_sutherland_data_ascii_default( PhysicalSet<SutherlandViscosity<NumericType>, ChemicalMixture<NumericType> >& mu )
  {
    antioch_deprecated();
    read_sutherland_data_ascii(mu, DefaultFilename::sutherland_data());
  }

} // end namespace Antioch

#endif // ANTIOCH_SUTHERLAND_PARSING_H
