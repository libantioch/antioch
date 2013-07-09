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
// $Id: eigen_utils.h 37170 2013-02-19 21:40:39Z roystgnr $
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef ANTIOCH_EIGEN_UTILS_DECL_H
#define ANTIOCH_EIGEN_UTILS_DECL_H

#ifdef ANTIOCH_METAPROGRAMMING_H
#  error eigen_utils_decl.h must be included before metaprogramming.h
#endif

#include "antioch_config.h"

#ifdef ANTIOCH_HAVE_EIGEN
// While this logic is Eigen-specific, these forward declarations deliberately
// do not require <Eigen/Dense> to be included.  This choice permits mixing
// header-only Eigen with header-only Antioch without configure-time flags.
#endif

// Forward declarations
namespace Eigen {
template<typename Derived>
class ArrayBase;

template<typename Derived>
class EigenBase;

template<typename Derived, typename ThenDerived, typename ElseDerived>
class Select;
}

#include "antioch/metaprogramming_decl.h"

namespace std
{

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
inline
_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
max(const _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& a,
    const _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& b);

} // end namespace std


namespace Antioch
{

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
inline
_Scalar
max(const _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& in);

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
struct has_size<_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> >;

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
struct size_type<_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> >;

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
struct value_type<_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> >;

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
  typename Scalar
>
inline
_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
constant_clone(const _Matrix<_Scalar, _Rows, _Cols, _Options,
	       _MaxRows, _MaxCols>& example,
	       const Scalar& value);

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
inline
void set_zero(_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& a);


template <typename Derived, typename ThenDerived, typename ElseDerived>
inline
const Eigen::Select<
  typename Eigen::ArrayBase<Derived>::PlainObject,
  typename Eigen::ArrayBase<ThenDerived>::PlainObject,
  typename Eigen::ArrayBase<ElseDerived>::PlainObject
>
if_else(
const Eigen::EigenBase<Derived>& condition,
const Eigen::EigenBase<ThenDerived>& if_true,
const Eigen::EigenBase<ElseDerived>& if_false);


} // end namespace Antioch

#endif //ANTIOCH_EIGEN_UTILS_DECL_H
