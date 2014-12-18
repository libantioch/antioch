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

#ifndef ANTIOCH_CEA_MIXTURE_ASCII_PARSING_H
#define ANTIOCH_CEA_MIXTURE_ASCII_PARSING_H

// C++
#include <iostream>
#include <vector>

// Antioch
#include "antioch/nasa_mixture_parsing.h"


/*! Everything here is deprecated, it is for backward compatibility,
    using the deprecated name/object CEAThermo...<NumericType>.

    We required to provide:
      - read_cea_mixture_data_ascii()
      - read_cea_mixture_data_ascii_default()

    for both descriptions (Mixture, dynamics).
*/

namespace Antioch
{
  // forward declarations
  template <typename NumericType>
  class CEACurveFit;

  template <typename NumericType, typename CurveType>
  class NASAThermoMixture;

  template <typename NumericType>
  class CEAThermodynamics;

  // required overload for backward compatibility
  template<class NumericType>
  void read_cea_mixture_data_ascii( CEAThermodynamics<NumericType>& thermo, const std::string &filename );

  template<class NumericType>
  void read_cea_mixture_data_ascii( NASAThermoMixture<NumericType,CEACurveFit<NumericType> >& thermo, const std::string &filename );

  template<class NumericType>
  void read_cea_mixture_data_ascii_default( CEAThermodynamics<NumericType>& thermo );

  template<class NumericType>
  void read_cea_mixture_data_ascii_default( NASAThermoMixture<NumericType, CEACurveFit<NumericType> >& thermo );

/* ------------------------ backward compatibility ---------------------*/


  template<class NumericType>
  void read_cea_mixture_data_ascii_default( CEAThermodynamics<NumericType>& thermo )
  {
    antioch_deprecated();
    read_cea_mixture_data_ascii(thermo, DefaultFilename::thermo_data());
  }

  template<class NumericType>
  void read_cea_mixture_data_ascii_default( NASAThermoMixture<NumericType, CEACurveFit<NumericType> >& thermo )
  {
    antioch_deprecated();
    read_cea_mixture_data_ascii(thermo, DefaultFilename::thermo_data());
  }

  template<class NumericType>
  void read_cea_mixture_data_ascii( CEAThermodynamics<NumericType> & thermo, const std::string &filename )
  {
    antioch_deprecated();
    ASCIIParser<NumericType> parser(filename);
    parser.read_thermodynamic_data(thermo);

   // Make sure we actually populated everything
    if( !thermo.check() )
    {
       std::cerr << "Error: CEA table not fully populated" << std::endl;
       antioch_error();
    }

    return;

  }

  template<class NumericType>
  void read_cea_mixture_data_ascii( NASAThermoMixture<NumericType, CEACurveFit<NumericType> > & thermo, const std::string &filename )
  {
    antioch_deprecated();
    ASCIIParser<NumericType> parser(filename);
    parser.read_thermodynamic_data(thermo);

   // Make sure we actually populated everything
    if( !thermo.check() )
    {
       std::cerr << "Error: CEA table not fully populated" << std::endl;
       antioch_error();
    }

    return;

  }


} // end namespace Antioch

#endif // ANTIOCH_CEA_MIXTURE_ASCII_PARSING_H
