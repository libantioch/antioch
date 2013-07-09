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
// $Id: metaprogramming.h 37170 2013-02-19 21:40:39Z roystgnr $
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef ANTIOCH_METAPROGRAMMING_DECL_H
#define ANTIOCH_METAPROGRAMMING_DECL_H

namespace Antioch
{
  // Helper metafunctions
  template <bool B, class T = void>
  struct enable_if_c {
    typedef T type;
  };

  template <class T>
  struct enable_if_c<false, T> {};

  // True for vector classes with size() methods
  template <typename T>
  struct has_size
  {
    const static bool value = false;
  };

  // size_type<T>::size is defined to be the result type of the size()
  // method for vector classes
  template <typename T>
  struct size_type { };


  // A class for uniformly assigning third-party vectorized numeric
  // types from scalar numeric types.  First-party vectorized numeric
  // types should correctly implement operator=(scalar)...
  template <typename T>
  struct value_type
  {
    typedef T type;

    typedef T raw_type;
  };

  template <typename T>
  inline
  T
  max (const T& in);

  template <typename T>
  struct value_type<const T>
  {
    typedef const typename value_type<T>::type type;

    typedef const typename value_type<T>::raw_type raw_type;
  };

  // A function for zero-initializing vectorized numeric types
  // while resizing them to match the example input
  template <typename T>
  inline
  T zero_clone(const T& example);

  // A function for initializing vectorized numeric types to a
  // constant // while resizing them to match the example input
  template <typename T, typename Scalar>
  inline
  T constant_clone(const T& example, const Scalar& value);

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

} // end namespace Antioch

#endif //ANTIOCH_METAPROGRAMMING_DECL_H
