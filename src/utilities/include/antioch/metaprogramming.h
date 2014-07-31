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

template <typename T, typename VectorScalar>
inline
T custom_clone(const T& /*example*/, const VectorScalar& values, unsigned int indexes) {return values[indexes];}

template <typename T, typename VectorScalar>
inline
T custom_clone(const T& /*example*/, const VectorScalar& values, const typename Antioch::rebind<T,unsigned int>::type & indexes)
{
  T returnval(indexes.size()); //bof bof - metaphysicl has size within type, this suppose that size_type<indexes>::type can be changed in value_type<T>::type
  
  for(std::size_t i = 0; i < indexes.size(); i++)
  {
      returnval[i] = values[indexes[i]];
  }
  return returnval;
}


template <typename VectorT>
inline
typename value_type<VectorT>::type
eval_index(const VectorT& vec, unsigned int index)
{
  return vec[index];
}

inline
bool conjunction(const bool & vec)
{
  return vec;
}

// FIXME - vexcl needs better 
// test: something in the flavor of
// vex::vector_expression<Condition>

template <typename T>
inline
bool conjunction(const T & vec)
{
  for(unsigned int i = 0; i < vec.size(); i++)
  {
    if(!vec[i])return false;
  }
  return true;
}

inline
bool disjunction(const bool & vec)
{
  return vec;
}

template <typename T>
inline
bool disjunction(const T & vec)
{
  for(unsigned int i = 0; i < vec.size(); i++)
  {
    if(vec[i])return true;
  }
  return false;
}

} // end namespace Antioch

#endif //ANTIOCH_METAPROGRAMMING_H
