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

#ifndef ANTIOCH_EIGEN_UTILS_H
#define ANTIOCH_EIGEN_UTILS_H

#ifdef ANTIOCH_METAPROGRAMMING_H
#  ifndef ANTIOCH_EIGEN_UTILS_DECL_H
#    error eigen_utils_decl.h must be included before metaprogramming.h
#  endif
#endif

#include "antioch_config.h"

#ifdef ANTIOCH_HAVE_EIGEN
// Though the following implementations are all valid without <Eigen/Dense>,
// successfully using them with Eigen types requires Eigen be included first.
// Configure-time Eigen support enforces this constraint but header-only Eigen
// may be mixed with header-only Antioch without configure-time flags.
#include <Eigen/Dense>
#endif

#include "antioch/metaprogramming.h"

// Notice _Matrix template templates might be Eigen::Matrix or Eigen::Array.
// Therefore, always use .array()- or .matrix()-like operations for robustness.
// Otherwise Eigen will complain with YOU_CANNOT_MIX_ARRAYS_AND_MATRICES.

namespace std
{

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
inline
_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
max(const _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& a,
    const _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& b)
{
  return a.array().max(b.array());
}

} // end namespace std


namespace Antioch
{

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
inline
_Scalar
max(const _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& in)
{
  return in.maxCoeff();
}

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
inline
_Scalar
min(const _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& in)
{
  return in.minCoeff();
}

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
struct has_size<_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> >
{
  static const bool value = true;
};

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
struct size_type<_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> >
{
  typedef typename
    _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>::Index
      type;
};

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
struct value_type<_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> >
{
  typedef _Matrix<
      _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols
    > container_type;
  typedef _Scalar type;
  typedef typename value_type<_Scalar>::raw_type raw_type;
};

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
inline
_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
zero_clone(const _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& ex)
{
  // We can't just use setZero here with arbitrary _Scalar types
  if (ex.size())
    return _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>(
        ex.rows(), ex.cols()).setConstant(zero_clone(ex[0]));

  return _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>(
      ex.rows(), ex.cols());
}

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols,
  typename Scalar2
>
inline
void
zero_clone(_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& output,
           const _Matrix<Scalar2, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& ex)
{
  // We can't just use setZero here with arbitrary _Scalar types
  output.resize(ex.rows(), ex.cols());
  output.setConstant(zero_clone(ex[0]));
}


template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols,
  typename Scalar
>
inline
_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
constant_clone(const _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& ex,
	       const Scalar& value)
{
  // We can't just use setZero here with arbitrary _Scalar types
  if (ex.size())
    return _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>(
        ex.rows(), ex.cols()).setConstant(value);

  return _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>(
      ex.rows(), ex.cols());
}


template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
inline
void set_zero(_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& a)
{
  // We can't just use setZero here with arbitrary _Scalar types
  if (a.size())
    a.setConstant (zero_clone(a[0]));
}


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
const Eigen::EigenBase<ElseDerived>& if_false)
{
  // Work around for current Eigen version bug in which we can't
  // call select() with expressions, only full objects.
  typename Eigen::ArrayBase<Derived>::PlainObject
    condition_array = condition;
  typename Eigen::ArrayBase<ThenDerived>::PlainObject
    if_true_array = if_true;
  typename Eigen::ArrayBase<ElseDerived>::PlainObject
    if_false_array = if_false;

  return condition_array.select(if_true_array, if_false_array);
}


} // end namespace Antioch

#endif //ANTIOCH_EIGEN_UTILS_H
