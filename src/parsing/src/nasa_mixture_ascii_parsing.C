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

#include "antioch/nasa_mixture_ascii_parsing.h"

// Antioch
#include "antioch/parsing_enum.h"
#include "antioch/nasa_mixture_parsing.h"

namespace Antioch
{

  template<class NumericType>
  void read_nasa_mixture_data_ascii( NASAThermoMixture<NumericType,NASA7CurveFit<NumericType> >& thermo,
                                     const std::string &filename )
  {
     antioch_deprecated();
     read_nasa_mixture_data( thermo, filename, CHEMKIN, true);
  }

  // Instantiate
  template void read_nasa_mixture_data_ascii<float>( NASAThermoMixture<float,NASA7CurveFit<float> >&, const std::string& );
  template void read_nasa_mixture_data_ascii<double>( NASAThermoMixture<double,NASA7CurveFit<double> >&, const std::string& );
  template void read_nasa_mixture_data_ascii<long double>( NASAThermoMixture<long double,NASA7CurveFit<long double> >&, const std::string& );

} // end namespace Antioch
