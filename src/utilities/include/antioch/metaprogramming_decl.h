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

#ifndef ANTIOCH_METAPROGRAMMING_DECL_H
#define ANTIOCH_METAPROGRAMMING_DECL_H

// C++
#include <cmath>   // For isnan

namespace Antioch
{
  // A workaround for trying to use commas in arguments to macros
#define ANTIOCH_COMMA ,

  // Allow us to use auto when sufficiently complete C++11 support is
  // available, while falling back on preselected types when C++11
  // is not available.
#ifdef ANTIOCH_HAVE_AUTO_THIS
#define ANTIOCH_AUTO(Type) auto
#define ANTIOCH_RETURNEXPR(Type, Expr) \
  -> typename Antioch::if_else_type<Antioch::return_auto<Type>::value, \
     decltype(Expr), \
     Type>::type
#define ANTIOCH_AUTOFUNC(Type, Expr) \
  -> typename Antioch::if_else_type<Antioch::return_auto<Type>::value, \
     decltype(Expr), \
     Type>::type \
  { return Expr; }
#else
#define ANTIOCH_AUTO(Type) Type
#define ANTIOCH_RETURNEXPR(Type, Expr)
#define ANTIOCH_AUTOFUNC(Type, Expr) { return Expr; }
#endif

  // Helper metafunctions
  template <bool B, class T = void>
  struct enable_if_c {
    typedef T type;
  };

  template <class T>
  struct enable_if_c<false, T> {};

  template <bool B, class T1, class T2>
  struct if_else_type {
    typedef T1 type;
  };

  template <class T1, class T2>
  struct if_else_type<false, T1, T2> {
    typedef T2 type;
  };

  // Base class to allow tag dispatching
  // http://www.generic-programming.org/languages/cpp/techniques.php
  struct numeric_library_tag {};

  // ::value == true for vector classes with size() methods
  template <typename T, typename Enable=void>
  struct has_size;

  // ::value == true for classes with expression templates that can be
  // safely returned via auto user functions
  template <typename T, typename Enable=void>
  struct return_auto;

  // size_type<T>::size is defined to be the result type of the size()
  // method for vector classes
  template <typename T, typename Enable=void>
  struct size_type;

  // A class to identify the underlying value type value_type<T>::type
  // of a numeric container type T.
  template <typename T, typename Enable=void>
  struct value_type;

  // A class to identify the concrete state type state_type<T>::type
  // of a (possibly template expression) type T.
  template <typename T, typename Enable=void>
  struct state_type;

  // A class to identify the underlying value type
  // raw_value_type<T>::type of
  // container-of-container types.
  template <typename T, typename Enable=void>
  struct raw_value_type;

  // A class for generating "int" from "float" or "valarray<int>" from
  // "valarray<float>" etc.
  template <typename Vector, typename NewScalar, typename Enable=void>
  struct rebind;

  /*
  template <typename T>
  inline
  typename value_type<T>::type
  max (const T& in);

  template <typename T>
  inline
  typename value_type<T>::type
  min (const T& in);
  */

  template <typename T, typename Enable=void>
  struct tag_type
  {
    typedef const numeric_library_tag type;
  };

  template <typename T>
  struct value_type<const T>
  {
    typedef const typename value_type<T>::type type;
  };

  template <typename T>
  struct raw_value_type<const T>
  {
    typedef const typename raw_value_type<T>::type type;
  };

  template <typename T1, typename T2>
  struct constructor_or_reference
  {
    typedef T1 type;
  };

  template <typename T>
  struct constructor_or_reference<T,T>
  {
    typedef T & type;
  };


  // A function for zero-initializing vectorized numeric types
  // while resizing them to match the example input
  template <typename T>
  inline
  T zero_clone(const T& example);

  // A function for zero-initializing vectorized numeric types
  // while resizing them to match the example input
  template <typename T1, typename T2>
  inline
  void zero_clone(T1& output, const T2& example);

  // A function for initializing vectorized numeric types to a
  // constant while resizing them to match the example input
  template <typename T, typename Scalar>
  inline
  T constant_clone(const T& example, const Scalar& value);

  // A function for initializing non vectorized numeric types to
  // custom constants stored in vector
  template <typename T, typename VectorScalar>
  inline
  T custom_clone(const T& example, const VectorScalar& values, unsigned int indexes);

  // A function for initializing vectorized numeric types to
  // custom constants stored in vector
  template <typename T, typename VectorScalar>
  inline
  T custom_clone(const T& example, const VectorScalar& values, const typename Antioch::rebind<T,unsigned int>::type & indexes);

  // A function for filling already-initialized vectorized numeric
  // types with a constant.
  template <typename T, typename Scalar>
  inline
  void constant_fill(T& output, const Scalar& value);

  // A function for initializing vectorized numeric types
  // while resizing them to match the example input
  template <typename T>
  inline
  void init_clone(T& output, const T& example);

  // A function for zero-setting vectorized numeric types
  template <typename T>
  inline
  void set_zero(T& output);

  // A function for initializing numeric vector types.  Resizes the
  // contents of vector-of-vector types but does not resize the outer
  // vector.
  template <typename Vector, typename Scalar>
  inline
  void init_constant(Vector& output, const Scalar& example);

  // A function for replicating the effect of the ternary operator on
  // numeric types.
  //
  // The first argument should be a boolean for a scalar type, or a
  // vector of booleans or an expression type for vector types.
  //
  // The second and third arguments should be expression types or
  // vector types.
  template <typename T>
  inline
  T if_else(bool condition,
	    T if_true,
	    T if_false);

  // A function for indexing a vector type with an (integral) numeric
  // type.
  //
  // The first argument should be a vector of (scalar or vectorized)
  // numeric types; the second argument should be a similarly scalar
  // or vectorized integer type with values valid as indices to the
  // first argument.  The output will be a similarly scalar or
  // vectorized numeric type.
  template <typename VectorT>
  inline
  typename value_type<VectorT>::type
  eval_index(const VectorT& vec, unsigned int index);

  //Root function for conjunction
  template <typename T>
  inline
  bool conjunction_root(const T & vec, numeric_library_tag);

  //Root function for disjunction
  template <typename T>
  inline
  bool disjunction_root(const T & vec, numeric_library_tag);

  // A function to obtain the conjunction of boolean
  template <typename T>
  inline
  bool conjunction(const T & vec);

  // A function to obtain the disjunction of boolean
  template <typename T>
  inline
  bool disjunction(const T & vec);

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
inline \
bool \
has_nan (const Type & in) \
{ \
  using std::isnan; \
 \
  return isnan(in); \
} \
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
ANTIOCH_PLAIN_SCALAR(bool);
ANTIOCH_PLAIN_SCALAR(short);
ANTIOCH_PLAIN_SCALAR(int);
ANTIOCH_PLAIN_SCALAR(long);
ANTIOCH_PLAIN_SCALAR(unsigned short);
ANTIOCH_PLAIN_SCALAR(unsigned int);
ANTIOCH_PLAIN_SCALAR(unsigned long);

} // end namespace Antioch

#endif //ANTIOCH_METAPROGRAMMING_DECL_H
