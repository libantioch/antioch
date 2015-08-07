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

#include "antioch/sutherland_parsing.h"

// Antioch
#include "antioch/sutherland_viscosity.h"
#include "antioch/mixture_viscosity.h"
#include "antioch/input_utils.h"
#include "antioch/sutherland_parsing_instantiate_macro.h"

// C++
#include <fstream>
#include <iostream>

namespace Antioch
{
  template<class NumericType>
  void read_sutherland_data_ascii( MixtureViscosity<SutherlandViscosity<NumericType>,NumericType >& mu,
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
	    const ChemicalMixture<NumericType>& chem_mixture = mu.chemical_mixture();

	    // Check if this is a species we want.
	    if( chem_mixture.species_name_map().find(name) !=
		chem_mixture.species_name_map().end() )
	      {
		// Pack up coefficients
		std::vector<NumericType> coeffs(2);
		coeffs[0] = a;
		coeffs[1] = b;
		mu.add(name, coeffs);
	      }
	  }
      }
    in.close();
  }

  template<class NumericType>
  void read_sutherland_data_ascii_default( MixtureViscosity<SutherlandViscosity<NumericType>,NumericType >& mu )
  {
    read_sutherland_data_ascii(mu, DefaultFilename::sutherland_data());
  }

  // Instantiate
  ANTIOCH_SUTHERLAND_PARSING_INSTANTIATE();

} // end namespace Antioch
