//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014 Paul T. Bauman, Sylvain Plessis, Roy H. Stonger
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
// $Id: metaprogramming.h 37170 2013-02-19 21:40:39Z roystgnr $
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef ANTIOCH_METAPROGRAMMING_H
#define ANTIOCH_METAPROGRAMMING_H

// Antioch
#include "antioch/metaprogramming_decl.h"

// C++
#include <cstddef> // For size_t

namespace Antioch
{
#define ANTIOCH_PLAIN_SCALAR(Type) \
 \
inline \
Type \
max (const Type& in) { return in; } \
 \
inline \
Type \
min (const Type& in) { return in; } \
 \
template <> \
inline \
Type if_else(bool condition, \
             Type if_true, \
             Type if_false) \
{ \
  if (condition) \
    return if_true; \
  else \
    return if_false; \
} \
 \
template <> \
struct has_size<Type> \
{ \
  const static bool value = false; \
}; \
 \
template <> \
struct return_auto<Type> \
{ \
  const static bool value = false; \
}; \
 \
template <> \
struct value_type<Type> \
{ \
  typedef Type type; \
}; \
 \
template <> \
struct raw_value_type<Type> \
{ \
  typedef Type type; \
}; \
 \
template <typename NewScalar> \
struct rebind<Type, NewScalar> \
{ \
  typedef NewScalar type; \
}


ANTIOCH_PLAIN_SCALAR(float);
ANTIOCH_PLAIN_SCALAR(double);
ANTIOCH_PLAIN_SCALAR(long double);
ANTIOCH_PLAIN_SCALAR(short);
ANTIOCH_PLAIN_SCALAR(int);
ANTIOCH_PLAIN_SCALAR(long);
ANTIOCH_PLAIN_SCALAR(unsigned short);
ANTIOCH_PLAIN_SCALAR(unsigned int);
ANTIOCH_PLAIN_SCALAR(unsigned long);

template <typename T>
inline
void set_zero(T& output) { output = 0; }

template <typename T>
inline
T zero_clone(const T& /* example */) { return 0; }

template <typename T, typename T2>
inline
void zero_clone(T& output, const T2& /* example */) { output = 0; }

template <typename T>
inline
void init_clone(T& output, const T& example) { output = example; }

template <typename Vector, typename Scalar>
inline
void init_constant(Vector& output, const Scalar& example)
{
  for (typename Antioch::size_type<Vector>::type i=0;
       i != output.size(); ++i)
    init_clone(output[i], example);
}

template <typename T, typename Scalar>
inline
T constant_clone(const T& /* example */, const Scalar& value) { return value; }

template <typename T, typename Scalar>
inline
void constant_fill(T& output, const Scalar& value) { output = value; }

} // end namespace Antioch

#endif //ANTIOCH_METAPROGRAMMING_H
