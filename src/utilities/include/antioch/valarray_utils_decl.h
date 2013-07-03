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

#ifndef ANTIOCH_VALARRAY_UTILS_DECL_H
#define ANTIOCH_VALARRAY_UTILS_DECL_H

#ifdef ANTIOCH_METAPROGRAMMING_H
#  error valarray_utils_decl.h must be included before metaprogramming.h
#endif

// Antioch
#include "antioch/metaprogramming_decl.h"

// C++
#include <iostream>
#include <valarray>

// Add some overloads that are blatantly missing from the std:: namespace

namespace std
{

template <typename T>
inline
std::ostream&
operator<< (std::ostream& output, const std::valarray<T>& a);

template <typename T, typename T2>
inline
std::valarray<T>
pow (const std::valarray<T>& in, const T2& n);

template <typename T>
inline
std::valarray<T>
max (const std::valarray<T>& a, const std::valarray<T>& b);

} // end namespace std


namespace Antioch
{

template <typename T>
inline
T
max (const std::valarray<T>& in);

template <typename T>
struct has_size<std::valarray<T> >;

template <typename T>
struct size_type<std::valarray<T> >;

template <typename T>
struct value_type<std::valarray<T> >;

template <typename T>
inline
std::valarray<T>
zero_clone(const std::valarray<T>& example);

template <typename T, typename Scalar>
inline
std::valarray<T>
constant_clone(const std::valarray<T>& example, const Scalar& value);

template <typename T>
inline
void
init_clone(std::valarray<T>& output, const std::valarray<T>& example);

} // end namespace Antioch


#endif //ANTIOCH_VALARRAY_UTILS_DECL_H
