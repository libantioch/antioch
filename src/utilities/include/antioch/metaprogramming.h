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

#ifndef ANTIOCH_METAPROGRAMMING_H
#define ANTIOCH_METAPROGRAMMING_H

// Antioch
#include "antioch/metaprogramming_decl.h"

// C++
#include <cstddef> // For size_t

namespace Antioch
{
  template <typename T>
  inline
  T
  max (const T& in)
  {
    return in;
  }

  template <typename T>
  inline
  T
  min (const T& in)
  {
    return in;
  }

  // A function for zero-initializing vectorized numeric types
  // while resizing them to match the example input
  template <typename T>
  inline
  T zero_clone(const T& /* example */) { return 0; }

  // A function for zero-initializing vectorized numeric types
  // while resizing them to match the example input
  template <typename T1, typename T2>
  inline
  void zero_clone(T1& output, const T2& /* example */) { output = 0; }

  // A function for initializing vectorized numeric types to a
  // constant // while resizing them to match the example input
  template <typename T, typename Scalar>
  inline
  T constant_clone(const T& example, const Scalar& value) { return value; }

  // A function for initializing vectorized numeric types
  // while resizing them to match the example input
  template <typename T>
  inline
  void init_clone(T& output, const T& example) { output = example; }

  // A function for zero-setting vectorized numeric types
  template <typename T>
  inline
  void set_zero(T& output) { output = 0; }

  // A function for initializing numeric vector types.  Resizes the
  // contents of vector-of-vector types but does not resize the outer
  // vector.
  template <typename Vector, typename Scalar>
  inline
  void init_constant(Vector& output, const Scalar& example)
  {
    // We can't just use setZero here with arbitrary Scalar types
    for (typename Antioch::size_type<Vector>::type i=0;
	 i != output.size(); ++i)
      init_clone(output[i], example);
  }

  // A default implementation for built-in types
  template <typename T>
  inline
  T if_else(bool condition,
	    T if_true,
	    T if_false)
  {
    if (condition)
      return if_true;
    else
      return if_false;
  }



} // end namespace Antioch

#endif //ANTIOCH_METAPROGRAMMING_H
