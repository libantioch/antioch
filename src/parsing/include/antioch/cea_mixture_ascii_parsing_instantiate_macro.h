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

#ifndef ANTIOCH_CEA_MIXTURE_ASCII_PARSING_INSTANTIATE_MACRO_H
#define ANTIOCH_CEA_MIXTURE_ASCII_PARSING_INSTANTIATE_MACRO_H

#define ANTIOCH_CEA_MIXTURE_DATA_ASCII_PARSING_TYPE_INSTANTIATE(type) \
  template void read_cea_mixture_data_ascii_default<type>( CEAThermoMixture<type>& ); \
  template void read_cea_mixture_data_ascii<type>( CEAThermoMixture<type>&, const std::string& ); \

#define ANTIOCH_CEA_MIXTURE_DATA_ASCII_PARSING_INSTANTIATE() \
  ANTIOCH_CEA_MIXTURE_DATA_ASCII_PARSING_TYPE_INSTANTIATE(float); \
  ANTIOCH_CEA_MIXTURE_DATA_ASCII_PARSING_TYPE_INSTANTIATE(double); \
  ANTIOCH_CEA_MIXTURE_DATA_ASCII_PARSING_TYPE_INSTANTIATE(long double)

#endif // ANTIOCH_CEA_MIXTURE_ASCII_PARSING_INSTANTIATE_MACRO_H
