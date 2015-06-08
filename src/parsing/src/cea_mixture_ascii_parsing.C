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

#include "antioch/cea_mixture_ascii_parsing.h"

// Antioch
#include "antioch/cea_mixture.h"
#include "antioch/ascii_parser.h"
#include "antioch/cea_mixture_ascii_parsing_instantiate_macro.h"

namespace Antioch
{

  template<class NumericType>
  void read_cea_mixture_data_ascii_default( CEAThermoMixture<NumericType >& thermo )
  {
    antioch_deprecated();
    read_cea_mixture_data(thermo, DefaultFilename::thermo_data(), ASCII, true);
  }

  template<class NumericType>
  void read_cea_mixture_data_ascii( CEAThermoMixture<NumericType > & thermo, const std::string &filename )
  {
    antioch_deprecated();
    read_cea_mixture_data( thermo, filename, ASCII, true );
  }

  template<class NumericType>
  void read_cea_mixture_data_ascii( CEAThermodynamics<NumericType > & thermo, const std::string &filename )
  {
    antioch_deprecated();
    ASCIIParser<NumericType> parser(filename);
    parser.read_thermodynamic_data(thermo);
  }

  ANTIOCH_CEA_MIXTURE_DATA_ASCII_PARSING_INSTANTIATE();

} // end namespace Antioch
