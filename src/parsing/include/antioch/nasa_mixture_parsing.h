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
#include <string>

// Antioch
#include "antioch/parsing_enum.h"
#include "antioch/default_filename.h"

namespace Antioch
{
  // Forward declarations
  template <class NumericType, class NASAFit>
  class NASAThermoMixture;

  // definition problems du to backward compatibilities...
  /*
     ASCIIParser -> CEAThermodynamics
                 ->  "antioch/cea_mixture_ascii_parsing.h"
                 ->  "antioch/cea_mixture_parsing.h"
                 ->  "antioch/nasa_mixture_parsing.h"
                 -> ASCIIParser
   */
  template <class NumericType>
  class ASCIIParser;

  template<class NumericType, typename CurveType >
  void read_nasa_mixture_data( NASAThermoMixture<NumericType, CurveType > & thermo, const std::string &filename = DefaultSourceFilename::thermo_data(), ParsingType = ASCII, bool verbose = true );

} // end namespace Antioch

#endif // ANTIOCH_NASA_MIXTURE_ASCII_PARSING_H
