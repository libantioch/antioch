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

#ifndef ANTIOCH_VEXCL_UTILS_H
#define ANTIOCH_VEXCL_UTILS_H

#ifdef ANTIOCH_METAPROGRAMMING_H
#  ifndef ANTIOCH_VEXCL_UTILS_DECL_H
#    error valarray_utils_decl.h must be included before metaprogramming.h
#  endif
#endif

// Antioch
#include "antioch/metaprogramming.h"
#include "antioch/antioch_asserts.h"

// C++
#include <cmath>
#include <iostream>
#include <valarray>


namespace Antioch
{

template <typename T>
inline
T
max (const vex::vector<T>& in)
{
  vex::Reductor<double, vex::MAX> max(in.queue_list());
  return max(in);
}

template <typename T>
struct has_size<vex::vector<T> >
{
  static const bool value = true;
};

template <typename T>
struct size_type<vex::vector<T> >
{
  typedef std::size_t type;
};

template <typename T>
struct value_type<vex::vector<T> >
{
  typedef vex::vector<T> container_type;
  typedef T type;
  typedef typename value_type<T>::raw_type raw_type;
};

template <typename T>
inline
vex::vector<T>
zero_clone(const vex::vector<T>& example)
{
  return vex::vector<T>(example.queue_list(), example.size());
}


template <typename T>
inline
void
init_clone(vex::vector<T>& output, const vex::vector<T>& example)
{
  // No resize method in VexCL!
  antioch_not_implemented();
}


} // end namespace Antioch


#endif //ANTIOCH_VEXCL_UTILS_H
