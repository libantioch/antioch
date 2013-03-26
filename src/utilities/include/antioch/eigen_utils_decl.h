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

#include "antioch_config.h"

#ifdef ANTIOCH_HAVE_EIGEN

#include "antioch/metaprogramming_decl.h"

#include <Eigen/Dense>

namespace std {

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline
Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
max (const Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& a,
     const Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& b);


template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline
Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
max (const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& a,
     const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& b);

}


namespace Antioch
{

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline
_Scalar
max (const Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& in);

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
struct value_type<Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> >;


template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline
Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
zero_clone(const Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& example);

// A function for zero-setting vectorized numeric types
template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline
void set_zero(Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& a);


template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline
_Scalar
max (const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& in);

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
struct value_type<Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> >;


template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline
Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
zero_clone(const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& example);

// A function for zero-setting vectorized numeric types
template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline
void set_zero(Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& a);


} // end namespace Antioch

#endif // ANTIOCH_HAVE_EIGEN

#endif //ANTIOCH_EIGEN_UTILS_DECL_H
