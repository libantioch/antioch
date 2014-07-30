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
ANTIOCH_AUTO(_Matrix<_Scalar ANTIOCH_COMMA  _Rows ANTIOCH_COMMA  _Cols ANTIOCH_COMMA  _Options ANTIOCH_COMMA  _MaxRows ANTIOCH_COMMA  _MaxCols>)
max(const _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & a,
    const _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & b)
ANTIOCH_AUTOFUNC(_Matrix<_Scalar ANTIOCH_COMMA  _Rows ANTIOCH_COMMA  _Cols ANTIOCH_COMMA  _Options ANTIOCH_COMMA  _MaxRows ANTIOCH_COMMA  _MaxCols>,
a.array().max(b.array()))

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
inline
ANTIOCH_AUTO(_Matrix<_Scalar ANTIOCH_COMMA  _Rows ANTIOCH_COMMA  _Cols ANTIOCH_COMMA  _Options ANTIOCH_COMMA  _MaxRows ANTIOCH_COMMA  _MaxCols>)
min(const _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & a,
    const _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & b)
ANTIOCH_AUTOFUNC(_Matrix<_Scalar ANTIOCH_COMMA  _Rows ANTIOCH_COMMA  _Cols ANTIOCH_COMMA  _Options ANTIOCH_COMMA  _MaxRows ANTIOCH_COMMA  _MaxCols>,
a.array().min(b.array()))

} // end namespace std


namespace Antioch
{

template <typename T>
inline
typename Antioch::enable_if_c<is_eigen<T>::value, 
  typename value_type<T>::type
  >::type
max(const T& in)
{
  return in.maxCoeff();
}

template <typename T>
inline
typename Antioch::enable_if_c<is_eigen<T>::value, 
  typename value_type<T>::type
  >::type
min(const T& in)
{
  return in.minCoeff();
}

template <typename T>
struct has_size<T, typename Antioch::enable_if_c<is_eigen<T>::value,void>::type>
{
  const static bool value = true;
};

template <typename T>
struct return_auto<T, typename Antioch::enable_if_c<is_eigen<T>::value,void>::type>
{
  const static bool value = true;
};

template <typename T>
struct size_type<T, typename Antioch::enable_if_c<is_eigen<T>::value,void>::type>
{
  typedef typename T::Index type;
};

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
struct value_type<_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> >
{
  typedef _Scalar type;
};

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
struct raw_value_type<_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> >
{
  typedef typename raw_value_type<_Scalar>::type type;
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
    return 
      _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
        (ex.rows(), ex.cols()).setConstant(zero_clone(ex[0]));

  return
    _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
      (ex.rows(), ex.cols());
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
constant_clone(const _Matrix<_Scalar, _Rows, _Cols, _Options,
	       _MaxRows, _MaxCols>& ex,
	       const Scalar& value)
{
  // We can't just use setZero here with arbitrary _Scalar types
  if (ex.size())
    return 
      _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
        (ex.rows(), ex.cols()).setConstant(value);

  return 
    _Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
      (ex.rows(), ex.cols());
}


template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols,
  typename Scalar
>
inline
void
constant_fill(_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows,
                      _MaxCols>& output,
	       const Scalar& value)
{
  output.fill(value);
}

template <
  template <typename, int, int, int, int, int> class _Matrix,
  typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
inline
void
set_zero(_Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& a)
{
  // We can't just use setZero here with arbitrary value_types
  if (a.size())
    a.setConstant (zero_clone(a[0]));
}

/*
template <template <typename, int, int, int, int, int> class _Matrix,
          typename T, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols,
  typename VectorScalar
>
inline
_Matrix<T,_Rows,_Cols,_Options,_MaxRows,_MaxCols> 
  custom_clone(const _Matrix<T,_Rows,_Cols,_Options,_MaxRows,_MaxCols> & /example/, const VectorScalar& values, const _Matrix<unsigned int,_Rows,_Cols,_Options,_MaxRows,_MaxCols> & indexes)
{
  _Matrix<T,_Rows,_Cols,_Options,_MaxRows,_MaxCols>  returnval;
  returnval.resize(indexes.size());
  for(std::size_t i = 0; i < indexes.size(); i++)
  {
      returnval[i] = values[indexes[i]];
  }
  return returnval;
}
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
const T2& if_false)
{
  return condition.select(if_true, if_false);
}


template <typename VectorT, 
  template <typename , int, int, int, int, int> class _Matrix,
  typename _UIntT, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols
>
inline
typename enable_if_c<
  is_eigen<typename value_type<VectorT>::type>::value,
  typename value_type<VectorT>::type
>::type
eval_index(const VectorT& vec, const _Matrix<_UIntT, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& index)
{
  typename value_type<VectorT>::type returnval = vec[0];
  for (std::size_t i=0; i != index.size(); ++i)
    returnval[i] = vec[index[i]][i];
  return returnval;
}


} // end namespace Antioch

#endif //ANTIOCH_EIGEN_UTILS_H
