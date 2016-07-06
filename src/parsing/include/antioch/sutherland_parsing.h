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


#ifndef ANTIOCH_SUTHERLAND_PARSING_H
#define ANTIOCH_SUTHERLAND_PARSING_H

// C++
#include <string>

namespace Antioch
{
  // Forward declarations
  template <typename NumericType>
  class SutherlandViscosity;

  template<typename Viscosity, typename NumericType>
  class MixtureViscosity;

  template<class NumericType>
  void read_sutherland_data_ascii( MixtureViscosity<SutherlandViscosity<NumericType>,NumericType >& mu,
				   const std::string &filename);

  template<class NumericType>
  void read_sutherland_data_ascii_default( MixtureViscosity<SutherlandViscosity<NumericType>,NumericType >& mu );

} // end namespace Antioch

#endif // ANTIOCH_SUTHERLAND_PARSING_H
