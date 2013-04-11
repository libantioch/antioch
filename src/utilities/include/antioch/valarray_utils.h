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

// Antioch
#include "antioch/metaprogramming.h"

// C++
#include <cmath>
#include <iostream>
#include <valarray>

// Add some overloads that are blatantly missing from the std:: namespace

namespace std
{

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

template <typename T, typename T2>
inline
valarray<T>
pow (const valarray<T>& in, const T2& n)
{
  valarray<T> out=in;
  const size_t size = in.size();
  for (size_t i=0; i != size; ++i)
    out[i] = pow(in[i], n);
  return out;
}

template <typename T>
inline
std::valarray<T>
max (const std::valarray<T>& a, const std::valarray<T>& b)
{
  using std::max;

  std::valarray<T> out=a;
  const size_t size = a.size();
  for (size_t i=0; i != size; ++i)
    out[i] = max(a[i], b[i]);
  return out;
}

} // end namespace std


namespace Antioch
{

template <typename T>
inline
T
max (const std::valarray<T>& in)
{
  return in.max();
}

template <typename T>
struct has_size<std::valarray<T> >
{
  static const bool value = true;
};

template <typename T>
struct value_type<std::valarray<T> >
{
  typedef std::valarray<T> container_type;
  typedef T type;
  typedef typename value_type<T>::raw_type raw_type;
};

template <typename T>
inline
std::valarray<T>
zero_clone(const std::valarray<T>& example)
{
  if (example.size())
    return std::valarray<T>(zero_clone(example[0]),example.size());
  else
    return std::valarray<T>();
}


template <typename T>
inline
void
init_clone(std::valarray<T>& output, const std::valarray<T>& example)
{
  const std::size_t sz = example.size();
  output.resize(sz);
  for (std::size_t i=0; i != sz; ++i)
    init_clone(output[i], example[i]);
}


} // end namespace Antioch


#endif //ANTIOCH_VALARRAY_UTILS_H
