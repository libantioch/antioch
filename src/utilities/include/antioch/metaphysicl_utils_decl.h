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

#ifndef ANTIOCH_METAPHYSICL_UTILS_DECL_H
#define ANTIOCH_METAPHYSICL_UTILS_DECL_H

#ifdef ANTIOCH_METAPROGRAMMING_H
#  error metaphysicl_utils_decl.h must be included before metaprogramming.h
#endif

#include "antioch_config.h"
#include "metaprogramming_decl.h"

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

#include "antioch/metaprogramming_decl.h"

// Specializations to match other Antioch workarounds

namespace Antioch
{

template <typename T>
struct is_metaphysicl {
  static const bool value = false;
};

template <std::size_t size, typename T>
struct is_metaphysicl<MetaPhysicL::NumberArray<size,T> > {
  static const bool value = true;
};

// Class to allow tag dispatching to MetaPhysicL specializations
struct metaphysicl_library_tag : public numeric_library_tag {};

// MetaPhysicL has no expression templates yet; all types store state
template <typename T>
struct state_type<T, typename enable_if_c<is_metaphysicl<T>::value,void>::type> {
  typedef T type;
};

template <std::size_t size, typename T, typename NewScalar>
struct rebind<MetaPhysicL::NumberArray<size,T>, NewScalar>
{
  typedef MetaPhysicL::NumberArray<size,NewScalar> type;
};

template <typename T>
inline
typename Antioch::enable_if_c<
  is_metaphysicl<T>::value,
  typename value_type<T>::type>::type
max (const T& in);

template <typename T>
inline
typename Antioch::enable_if_c<
  is_metaphysicl<T>::value,
  typename value_type<T>::type>::type
min (const T& in);

template <typename T>
inline
typename Antioch::enable_if_c<
  is_metaphysicl<T>::value,
  bool>::type
has_nan (const T& in);

template <typename T>
struct has_size<T, typename Antioch::enable_if_c<is_metaphysicl<T>::value,void>::type>;

template <typename T>
struct return_auto<T, typename Antioch::enable_if_c<is_metaphysicl<T>::value,void>::type>;

template <typename T>
struct size_type<T, typename Antioch::enable_if_c<is_metaphysicl<T>::value,void>::type>;

template <typename T>
struct value_type<T, typename Antioch::enable_if_c<is_metaphysicl<T>::value,void>::type>;

template <typename T>
struct raw_value_type<T, typename Antioch::enable_if_c<is_metaphysicl<T>::value,void>::type>;

template <typename Tbool, typename Ttrue, typename Tfalse>
inline
typename Antioch::enable_if_c<
  is_metaphysicl<Tbool>::value &&
  is_metaphysicl<Ttrue>::value &&
  is_metaphysicl<Tfalse>::value,
  typename state_type<Ttrue>::type>::type
if_else(const Tbool& condition,
        const Ttrue& if_true,
        const Tfalse& if_false);

template <typename VectorT, typename UIntType>
inline
typename Antioch::enable_if_c<
   is_metaphysicl<typename Antioch::value_type<VectorT>::type>::value &&
   is_metaphysicl<UIntType>::value,
   typename Antioch::value_type<VectorT>::type
>::type
eval_index(const VectorT & vec, const UIntType & indexes);
} // end namespace Antioch

#endif // ANTIOCH_METAPHYSICL_UTILS_DECL_H
