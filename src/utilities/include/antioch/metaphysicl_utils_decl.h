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

#ifndef ANTIOCH_METAPHYSICL_UTILS_DECL_H
#define ANTIOCH_METAPHYSICL_UTILS_DECL_H

#include "antioch_config.h"

#ifdef ANTIOCH_HAVE_METAPHYSICL
#include "metaphysicl/numberarray.h"

#include "antioch/metaprogramming_decl.h"

// Specializations to match other Antioch workarounds

namespace Antioch
{

template <std::size_t size, typename T>
inline
T
max (const MetaPhysicL::NumberArray<size,T>& in);

template <std::size_t size, typename T>
struct has_size<MetaPhysicL::NumberArray<size,T> >;

template <std::size_t size, typename T>
struct value_type<MetaPhysicL::NumberArray<size,T> >;

template <std::size_t size, typename T>
inline
MetaPhysicL::NumberArray<size,T>
zero_clone(const MetaPhysicL::NumberArray<size,T>& example);

} // end namespace Antioch

#endif // ANTIOCH_HAVE_METAPHYSICL

#endif // ANTIOCH_METAPHYSICL_UTILS_DECL_H
