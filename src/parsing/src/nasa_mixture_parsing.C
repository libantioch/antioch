//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014-2016 Paul T. Bauman, Benjamin S. Kirk,
//                         Sylvain Plessis, Roy H. Stonger
//
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

#include "antioch/nasa_mixture_parsing.h"
#include "antioch/nasa_mixture_ascii_parsing.h"

// Antioch
#include "antioch/nasa_mixture_parsing_instantiate_macro.h"
#include "antioch/parsing_enum.h"
#include "antioch/ascii_parser.h"
#include "antioch/xml_parser.h"
#include "antioch/chemkin_parser.h"
#include "antioch/nasa_mixture.h"

namespace Antioch
{
  template<class NumericType, typename CurveType>
  void read_nasa_mixture_data( NASAThermoMixture<NumericType, CurveType >& thermo, const std::string &filename, ParsingType type, bool verbose )
  {

    ParserBase<NumericType> * parser(NULL);
    switch(type)
      {
      case ASCII:
        parser = new ASCIIParser<NumericType>(filename,verbose);
        break;
      case CHEMKIN:
        parser = new ChemKinParser<NumericType>(filename,verbose);
        break;
      case XML:
        parser = new XMLParser<NumericType>(filename,verbose);
        break;
      default:
        antioch_parsing_error("unknown type");
      }

    parser->read_thermodynamic_data(thermo);

    // Make sure we actually populated everything
    if( !thermo.check() )
      {
	std::cerr << "Error: NASA table not fully populated" << std::endl;
	antioch_error();
      }

    return;
  }

  template<class NumericType>
  void read_nasa_mixture_data_ascii( NASAThermoMixture<NumericType, NASA7CurveFit<NumericType> >& thermo, const std::string &filename )
  {
    antioch_deprecated();
    read_nasa_mixture_data( thermo, filename, CHEMKIN, true);
  }

  // Instantiate
  ANTIOCH_NASA_MIXTURE_PARSING_INSTANTIATE(NASA7CurveFit);
  ANTIOCH_NASA_MIXTURE_PARSING_INSTANTIATE(NASA9CurveFit);
  ANTIOCH_NASA_MIXTURE_PARSING_INSTANTIATE(CEACurveFit);

  template void read_nasa_mixture_data_ascii<float>( NASAThermoMixture<float,NASA7CurveFit<float> >&, const std::string& );
  template void read_nasa_mixture_data_ascii<double>( NASAThermoMixture<double,NASA7CurveFit<double> >&, const std::string& );
  template void read_nasa_mixture_data_ascii<long double>( NASAThermoMixture<long double,NASA7CurveFit<long double> >&, const std::string& );

} // end namespace Antioch
