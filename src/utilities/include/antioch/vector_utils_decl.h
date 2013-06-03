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

#ifndef ANTIOCH_VECTOR_UTILS_DECL_H
#define ANTIOCH_VECTOR_UTILS_DECL_H

#ifdef ANTIOCH_METAPROGRAMMING_H
#  error vector_utils_decl.h must be included before metaprogramming.h
#endif

// Antioch
#include "antioch/metaprogramming_decl.h"

// C++
#include <iostream>
#include <vector>

// Add some overloads

namespace std
{

template <typename T>
inline
std::ostream&
operator<< (std::ostream& output, const std::vector<T>& a);

} // end namespace std

namespace Antioch
{

template <typename T>
struct has_size<std::vector<T> >;

template <typename T>
struct size_type<std::vector<T> >;

template <typename T>
struct value_type<std::vector<T> >;

template <typename T>
inline
std::vector<T>
zero_clone(const std::vector<T>& example);

// A function for zero-setting vectorized numeric types
template <typename T>
inline
void set_zero(std::vector<T>& a);

} // end namespace Antioch


#endif //ANTIOCH_VECTOR_UTILS_DECL_H
