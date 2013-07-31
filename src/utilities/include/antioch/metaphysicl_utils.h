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

#ifndef ANTIOCH_METAPHYSICL_UTILS_H
#define ANTIOCH_METAPHYSICL_UTILS_H

#ifdef ANTIOCH_METAPROGRAMMING_H
#  ifndef ANTIOCH_METAPHYSICL_UTILS_DECL_H
#    error metaphysicl_utils_decl.h must be included before metaprogramming.h
#  endif
#endif

#include "antioch_config.h"

#include <cstddef> // std::size_t

#ifdef ANTIOCH_HAVE_METAPHYSICL
// Though the following implementations are all valid without
// <metaphysicl/numberarray.h>, successfully using them with
// MetaPhysicL types requires MetaPhysicL to be included first.
// Configure-time MetaPhysicL support enforces this constraint but
// header-only MetaPhysicL may be mixed with header-only Antioch
// without configure-time flags.
#include "metaphysicl/numberarray.h"
#else
// Forward declaration instead
namespace MetaPhysicL {
template <std::size_t size, typename T> class NumberArray;
}
#endif

#include "antioch/metaprogramming.h"

// Specializations to match other Antioch workarounds

namespace Antioch
{

template <typename T>
inline
typename Antioch::enable_if_c<
  is_metaphysicl<T>::value,
  typename value_type<T>::type>::type
max (const T& in)
{
  using std::max;

  typename value_type<T>::type maxval = in[0];
  const std::size_t size = in.size();
  for (std::size_t i = 1; i < size; ++i)
    maxval = max(maxval, in[i]);

  return maxval;
}

template <typename T>
inline
typename Antioch::enable_if_c<
  is_metaphysicl<T>::value,
  typename value_type<T>::type>::type
min (const T& in)
{
  using std::min;

  typename value_type<T>::type minval = in[0];
  const std::size_t size = in.size();
  for (std::size_t i = 1; i < size; ++i)
    minval = min(minval, in[i]);

  return minval;
}

template <typename T>
struct has_size<T, typename Antioch::enable_if_c<is_metaphysicl<T>::value,void>::type>
{
  static const bool value = true;
};

template <typename T>
struct return_auto<T, typename Antioch::enable_if_c<is_metaphysicl<T>::value,void>::type>
{
  static const bool value = true;
};

template <typename T>
struct size_type<T, typename Antioch::enable_if_c<is_metaphysicl<T>::value,void>::type>
{
  typedef std::size_t type;
};

template <typename T>
struct value_type<T, typename Antioch::enable_if_c<is_metaphysicl<T>::value,void>::type>
{
  typedef typename T::value_type type;
};

template <typename T>
struct raw_value_type<T, typename Antioch::enable_if_c<is_metaphysicl<T>::value,void>::type>
{
  typedef typename raw_value_type<T>::type type;
};

template <typename Tbool, typename Ttrue, typename Tfalse>
inline
typename Antioch::enable_if_c<
  is_metaphysicl<Tbool>::value &&
  is_metaphysicl<Ttrue>::value &&
  is_metaphysicl<Tfalse>::value,
  typename state_type<Ttrue>::type>::type
if_else(const Tbool& condition,
        const Ttrue& if_true,
        const Tfalse& if_false)
{
  typename state_type<Ttrue>::type returnval;
  antioch_assert_equal_to(condition.size(), if_true.size());
  antioch_assert_equal_to(condition.size(), if_false.size());

  for (std::size_t i=0; i != condition.size(); ++i)
    returnval[i] = condition[i] ? if_true[i] : if_false[i];

  return returnval;
}

} // end namespace Antioch

#endif // ANTIOCH_METAPHYSICL_UTILS_H
