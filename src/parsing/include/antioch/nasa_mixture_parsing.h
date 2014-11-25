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
  class NASACurveFit;

  // New declarations

  template<class NumericType,typename Parser>
  void read_nasa_mixture_data( NASAThermoMixture<NumericType,NASACurveFit<NumericType> >& thermo, const std::string &filename );

 
  /* ------------------------- Inline Functions -------------------------*/
  template<class NumericType, typename Parser>
  inline
  void read_nasa_mixture_data( NASAThermoMixture<NumericType, NASACurveFit<NumericType> >& thermo, const std::string &filename )
  {
    
   Parser parser(filename);

   parser.read_thermodynamic_data(thermo);

   // Make sure we actually populated everything
   if( !thermo.check() )
     {
	std::cerr << "Error: NASA table not fully populated" << std::endl;
	antioch_error();
      }

    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_NASA_MIXTURE_ASCII_PARSING_H
