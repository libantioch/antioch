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

#ifndef ANTIOCH_XML_PARSER_INSTANTIATION_MACRO_H
#define ANTIOCH_XML_PARSER_INSTANTIATION_MACRO_H

#define ANTIOCH_XML_PARSER_TYPE_INSTANTIATE(type)                     \
  template void XMLParser<type>::read_thermodynamic_data_root<NASAThermoMixture<type,NASA7CurveFit<type> > >(NASAThermoMixture<type,NASA7CurveFit<type> >& ); \
  template void XMLParser<type>::read_thermodynamic_data_root<NASAThermoMixture<type,NASA9CurveFit<type> > >(NASAThermoMixture<type,NASA9CurveFit<type> >& ); \
  template void XMLParser<type>::read_thermodynamic_data_root<NASAThermoMixture<type,CEACurveFit<type> > >(NASAThermoMixture<type,CEACurveFit<type> >& ); \

#define ANTIOCH_XML_PARSER_INSTANTIATE() \
  ANTIOCH_XML_PARSER_TYPE_INSTANTIATE(float); \
  ANTIOCH_XML_PARSER_TYPE_INSTANTIATE(double); \
  ANTIOCH_XML_PARSER_TYPE_INSTANTIATE(long double)

#endif // ANTIOCH_XML_PARSER_INSTANTIATION_MACRO_H
