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

#ifdef ANTIOCH_METAPROGRAMMING_H
#  ifndef ANTIOCH_VALARRAY_UTILS_DECL_H
#    error valarray_utils_decl.h must be included before metaprogramming.h
#  endif
#endif

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
typename Antioch::enable_if_c<
  Antioch::is_valarray<T>::value,
  std::ostream&>::type
operator<< (std::ostream& output, const T& a)
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
typename Antioch::enable_if_c<
  Antioch::is_valarray<T>::value,
  typename Antioch::state_type<T>::type>::type
pow (const T& in, const T2& n)
{
  typename Antioch::state_type<T>::type out=in;
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

  const size_t size = a.size();
  std::valarray<T> out(size);
  for (size_t i=0; i != size; ++i)
    out[i] = max(a[i], b[i]);
  return out;
}


template <typename T>
inline
std::valarray<T>
min (const std::valarray<T>& a, const std::valarray<T>& b)
{
  using std::min;

  const size_t size = a.size();
  std::valarray<T> out(size);
  for (size_t i=0; i != size; ++i)
    out[i] = min(a[i], b[i]);
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
inline
T
min (const std::valarray<T>& in)
{
  return in.min();
}

template <typename T>
struct return_auto<T, typename Antioch::enable_if_c<is_valarray<T>::value,void>::type>
{
  static const bool value = false;
};

template <typename T>
struct has_size<std::valarray<T> >
{
  static const bool value = true;
};

template <typename T>
struct size_type<std::valarray<T> >
{
  typedef std::size_t type;
};

template <typename T>
struct value_type<std::valarray<T> >
{
  typedef T type;
};

template <typename T>
struct raw_value_type<std::valarray<T> >
{
  typedef typename raw_value_type<T>::type type;
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

template <typename T1, typename T2>
inline
void
zero_clone(std::valarray<T1>& output, const std::valarray<T2>& example)
{
  const std::size_t sz = example.size();
  output.resize(sz);
  for (std::size_t i=0; i != sz; ++i)
    Antioch::zero_clone(output[i], example[i]);
}

template <typename T, typename Scalar>
inline
std::valarray<T>
constant_clone(const std::valarray<T>& example, const Scalar& value)
{
  if (example.size())
    return std::valarray<T>(value,example.size());
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


template <typename T>
inline
std::valarray<T>
if_else(const std::valarray<bool>& condition,
        const std::valarray<T>& if_true,
        const std::valarray<T>& if_false)
{
  antioch_assert_equal_to(condition.size(), if_true.size());
  antioch_assert_equal_to(condition.size(), if_false.size());

  const std::size_t size = condition.size();
  std::valarray<T> returnval(size);

  for (std::size_t i=0; i != size; ++i)
    returnval[i] = condition[i] ? if_true[i] : if_false[i];

  return returnval;
}

template <typename VectorT, typename IntT>
inline
typename value_type<VectorT>::type
eval_index(const VectorT& vec, const std::valarray<IntT>& index)
{
  typename value_type<VectorT>::type returnval;
  std::size_t sz = index.size();
  returnval.resize(sz);
  for (std::size_t i=0; i != sz; ++i)
    returnval[i] = vec[index[i]][i];
  return returnval;
}

} // end namespace Antioch


#endif //ANTIOCH_VALARRAY_UTILS_H
