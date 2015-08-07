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

#include "antioch/cea_mixture_parsing.h"

// Antioch
#include "antioch/ascii_parser.h"
#include "antioch/chemkin_parser.h"
#include "antioch/xml_parser.h"
#include "antioch/cea_mixture.h"

namespace Antioch
{
  template<class NumericType>
  void read_cea_mixture_data( CEAThermoMixture<NumericType >& thermo, const std::string & filename, ParsingType type, bool verbose )
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
       std::cerr << "Error: CEA table not fully populated" << std::endl;
       antioch_error();
    }
  }

  // Instantiate
  template void read_cea_mixture_data<float>( CEAThermoMixture<float>&,const std::string&, ParsingType, bool );
  template void read_cea_mixture_data<double>( CEAThermoMixture<double>&,const std::string&, ParsingType, bool );
  template void read_cea_mixture_data<long double>( CEAThermoMixture<long double>&,const std::string&, ParsingType, bool );

} // end namespace Antioch
