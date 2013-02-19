//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// Antioch - A Gas Dynamics Thermochemistry Library
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
//
// $Id: valarray_utils.h 37170 2013-02-19 21:40:39Z roystgnr $
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef ANTIOCH_VALARRAY_UTILS_H
#define ANTIOCH_VALARRAY_UTILS_H

#include <cmath>
#include <iostream>
#include <valarray>

// Add some overloads that are blatantly missing from the std:: namespace

template <typename T>
inline
std::ostream&
operator<< (std::ostream& output, const std::valarray<T>& a)
{
  output << '{';
  const std::size_t size = a.size();
  if (size)
    output << a[0];
  for (std::size_t i=1; i<size; ++i)
    output << ',' << a[i];
  output << '}';
  return output;
}

namespace std
{

template <typename T>
inline
valarray<T>
pow (valarray<T> in, int n)
{
  const size_t size = in.size();
  for (unsigned int i=0; i != size; ++i)
    in[i] = pow(in[i], n);
  return in;
}

} // end namespace std

#endif //ANTIOCH_VALARRAY_UTILS_H
