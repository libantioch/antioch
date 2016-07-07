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

#ifndef ANTIOCH_EIGEN_UTILS_DECL_H
#define ANTIOCH_EIGEN_UTILS_DECL_H

#ifdef ANTIOCH_METAPROGRAMMING_H
#  error eigen_utils_decl.h must be included before metaprogramming.h
#endif

#include "antioch_config.h"
#include "metaprogramming_decl.h"

#include <type_traits> // std::enable_if

#ifdef ANTIOCH_HAVE_EIGEN
// While this logic is Eigen-specific, these forward declarations deliberately
// do not require <Eigen/Dense> to be included.  This choice permits mixing
// header-only Eigen with header-only Antioch without configure-time flags.
#endif

// Forward declarations
namespace Eigen {
template<typename Derived, typename ThenDerived, typename ElseDerived>
class Select;
}

// We want this code to apply to Eigen matrices, arrays, and
// expression templates

namespace Antioch {

template <typename T>
struct is_eigen {
  static const bool value = false;
};

// FIXME: need to add specializations for expression template types
template <typename T>
struct state_type<T, typename enable_if_c<is_eigen<T>::value,void>::type> {
  typedef T type;
};

// Class to allow tag dispatching to Eigen specializations
struct eigen_library_tag : public numeric_library_tag {};

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
struct is_eigen
<_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> >
{
  static const bool value = true;
};


}

#include "antioch/metaprogramming_decl.h"

namespace std
{

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
inline
ANTIOCH_AUTO(_Matrix<_Scalar ANTIOCH_COMMA  _Rows ANTIOCH_COMMA  _Cols ANTIOCH_COMMA  _Options ANTIOCH_COMMA  _MaxRows ANTIOCH_COMMA  _MaxCols>)
max(const _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & a,
    const _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & b)
ANTIOCH_RETURNEXPR(_Matrix<_Scalar ANTIOCH_COMMA  _Rows ANTIOCH_COMMA  _Cols ANTIOCH_COMMA  _Options ANTIOCH_COMMA  _MaxRows ANTIOCH_COMMA  _MaxCols>,
a.array().max(b.array()));

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
inline
ANTIOCH_AUTO(_Matrix<_Scalar ANTIOCH_COMMA  _Rows ANTIOCH_COMMA  _Cols ANTIOCH_COMMA  _Options ANTIOCH_COMMA  _MaxRows ANTIOCH_COMMA  _MaxCols>)
min(const _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & a,
    const _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & b)
ANTIOCH_RETURNEXPR(_Matrix<_Scalar ANTIOCH_COMMA  _Rows ANTIOCH_COMMA  _Cols ANTIOCH_COMMA  _Options ANTIOCH_COMMA  _MaxRows ANTIOCH_COMMA  _MaxCols>,
a.array().min(b.array()));

} // end namespace std


namespace Antioch
{
template <typename T>
struct tag_type<T,
                typename std::enable_if<is_eigen<T>::value>::type
>
{
  typedef eigen_library_tag type;
};

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols,
  typename NewScalar
>
struct rebind<_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>, NewScalar>
{
  typedef _Matrix<NewScalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> type;
};

template <typename T>
struct has_size<T, typename Antioch::enable_if_c<is_eigen<T>::value,void>::type>;

template <typename T>
struct return_auto<T, typename Antioch::enable_if_c<is_eigen<T>::value,void>::type>;

template <typename T>
struct size_type<T, typename Antioch::enable_if_c<is_eigen<T>::value,void>::type>;

// FIXME - we need value_type defined for expression templates too
template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
struct value_type<_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> >;

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
struct raw_value_type<_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> >;

template <typename T>
inline
typename Antioch::enable_if_c<is_eigen<T>::value,
  typename value_type<T>::type
  >::type
max(const T& in);

template <typename T>
inline
typename Antioch::enable_if_c<is_eigen<T>::value,
  typename value_type<T>::type
  >::type
min(const T& in);

template <typename T>
inline
typename Antioch::enable_if_c<is_eigen<T>::value,
  bool>::type
has_nan(const T& in);


template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
inline
_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
zero_clone(const _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& ex);

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols,
  typename Scalar2
>
inline
void
zero_clone(_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& output,
	   const _Matrix<Scalar2, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& ex);

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols,
  typename Scalar
>
inline
_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
constant_clone(const _Matrix<_Scalar, _Rows, _Cols, _Options,
	       _MaxRows, _MaxCols>& ex,
	       const Scalar& value);

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols,
  typename Scalar
>
inline
void
constant_fill(_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows,
                      _MaxCols>& output,
	       const Scalar& value);


template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
inline
void
set_zero(_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& a);

/*
template <template <typename, int, int, int, int, int> class _Matrix,
          typename T, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols,
  typename VectorScalar
>
inline
_Matrix<T,_Rows,_Cols,_Options,_MaxRows,_MaxCols>
  custom_clone(const _Matrix<T,_Rows,_Cols,_Options,_MaxRows,_MaxCols> & example, const VectorScalar& values, const _Matrix<unsigned int,_Rows,_Cols,_Options,_MaxRows,_MaxCols> & indexes);
*/

template <typename Condition, typename T1, typename T2>
inline
typename enable_if_c<
  is_eigen<T1>::value &&
  is_eigen<T2>::value,
  typename state_type<T1>::type
>::type
if_else(
const Condition& condition,
const T1& if_true,
const T2& if_false);


template <typename VectorT,
  template <typename, int, int, int, int, int> class _Matrix,
  typename _UIntT, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
inline
typename enable_if_c<
  is_eigen<typename value_type<VectorT>::type>::value,
  typename value_type<VectorT>::type
>::type
eval_index(const VectorT& vec, const _Matrix<_UIntT, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& index);

template <typename T>
inline
bool conjunction_root(const T & vec, eigen_library_tag);

template <typename T>
inline
bool disjunction_root(const T & vec, eigen_library_tag);

} // end namespace Antioch

#endif //ANTIOCH_EIGEN_UTILS_DECL_H
