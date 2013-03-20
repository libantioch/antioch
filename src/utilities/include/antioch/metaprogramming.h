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

namespace Antioch
{
  // Helper metafunctions
  template <bool B, class T = void>
  struct enable_if_c {
    typedef T type;
  };

  template <class T>
  struct enable_if_c<false, T> {};

  template <typename T>
  class has_size
  {
    typedef char no;
    typedef char yes[2];
    template <class C> static yes& test(char (*)[sizeof(&C::size)]);
    template <class C> static no& test(...);
  public:
    const static bool value = (sizeof(test<T>(0)) == sizeof(yes&));
  };

  // A class for uniformly assigning third-party vectorized numeric
  // types from scalar numeric types.  First-party vectorized numeric
  // types should correctly implement operator=(scalar)...
  template <typename T>
  struct value_type
  {
    typedef T type;

    static inline
    T constant(const type& in) { return in; }
  };

  template <typename T>
  inline
  T
  max (const T& in)
  {
    return in;
  }

  template <typename T>
  struct value_type<const T>
  {
    typedef const typename value_type<T>::type type;

    static inline
    T constant(const type& in) { return in; }
  };

  // A function for zero-initializing vectorized numeric types
  template <typename T>
  inline
  T zero_clone(const T& example) { return 0; }

  // A function for zero-setting vectorized numeric types
  template <typename T>
  inline
  void set_zero(T& a) { a = 0; }
} // end namespace Antioch

#endif //ANTIOCH_METAPROGRAMMING_H
