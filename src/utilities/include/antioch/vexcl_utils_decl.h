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

#ifndef ANTIOCH_VEXCL_UTILS_DECL_H
#define ANTIOCH_VEXCL_UTILS_DECL_H

#ifdef ANTIOCH_METAPROGRAMMING_H
#  error valarray_utils_decl.h must be included before metaprogramming.h
#endif

// Antioch
#include "antioch/metaprogramming_decl.h"

// C++
#include <iostream>
#include <valarray>

namespace Antioch
{

template <typename T>
inline
T
max (const vex::vector<T>& in);

template <typename T>
struct has_size<vex::vector<T> >;

template <typename T>
struct size_type<vex::vector<T> >;

template <typename T>
struct value_type<vex::vector<T> >;

template <typename T>
inline
vex::vector<T>
zero_clone(const vex::vector<T>& example);

template <typename T>
inline
void
init_clone(vex::vector<T>& output, const vex::vector<T>& example);

} // end namespace Antioch


#endif //ANTIOCH_VEXCL_UTILS_DECL_H
