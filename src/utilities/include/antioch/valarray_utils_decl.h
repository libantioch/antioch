//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014 Paul T. Bauman, Benjamin S. Kirk, Sylvain Plessis,
//                    Roy H. Stonger
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

namespace Antioch
{

template <typename T>
struct is_valarray {
  static const bool value = false;
};

template <typename T>
struct is_valarray<std::valarray<T> > {
  static const bool value = true;
};

// FIXME: need to add specializations for expression template types
template <typename T>
struct state_type<T, typename enable_if_c<is_valarray<T>::value,void>::type> {
  typedef T type;
};
// Class to allow tag dispatching to std::valarray specializations
struct valarray_library_tag : public numeric_library_tag {};

}


namespace std
{

template <typename T>
inline
typename Antioch::enable_if_c<
  Antioch::is_valarray<T>::value,
  std::ostream&>::type
operator<< (std::ostream& output, const T& a);

template <typename T, typename T2>
inline
typename Antioch::enable_if_c<
  Antioch::is_valarray<T>::value,
  typename Antioch::state_type<T>::type>::type
pow (const T& in, const T2& n);

template <typename T>
inline
std::valarray<T>
max (const std::valarray<T>& a, const std::valarray<T>& b);

template <typename T>
inline
std::valarray<T>
min (const std::valarray<T>& a, const std::valarray<T>& b);

} // end namespace std


namespace Antioch
{
template <typename T, typename NewScalar>
struct rebind<std::valarray<T>, NewScalar>
{
  typedef std::valarray<NewScalar> type;
};

template <typename T>
inline
T
max (const std::valarray<T>& in);

template <typename T>
inline
T
min (const std::valarray<T>& in);

template <typename T>
struct return_auto<T, typename Antioch::enable_if_c<is_valarray<T>::value,void>::type>;

template <typename T>
struct has_size<std::valarray<T> >;

template <typename T>
struct size_type<std::valarray<T> >;

template <typename T>
struct value_type<std::valarray<T> >;

template <typename T>
struct raw_value_type<std::valarray<T> >;

template <typename T>
inline
std::valarray<T>
zero_clone(const std::valarray<T>& example);

template <typename T1, typename T2>
inline
void
zero_clone(std::valarray<T1>& output, const std::valarray<T2>& example);

template <typename T, typename Scalar>
inline
std::valarray<T>
constant_clone(const std::valarray<T>& example, const Scalar& value);

template <typename T>
inline
void
init_clone(std::valarray<T>& output, const std::valarray<T>& example);

/*template <typename T, typename VectorScalar>
inline
std::valarray<T> custom_clone(const std::valarray<T>& example, const VectorScalar& values, const std::valarray<unsigned int> & indexes);
*/

template <typename T>
inline
std::valarray<T>
if_else(const std::valarray<bool>& condition,
        const std::valarray<T>& if_true,
        const std::valarray<T>& if_false);


template <typename VectorT>
inline
typename Antioch::enable_if_c<
        Antioch::is_valarray<typename value_type<VectorT>::type>::value,
        typename value_type<VectorT>::type
>::type
eval_index(const VectorT& vec, const std::valarray<unsigned int>& index);

} // end namespace Antioch


#endif //ANTIOCH_VALARRAY_UTILS_DECL_H
