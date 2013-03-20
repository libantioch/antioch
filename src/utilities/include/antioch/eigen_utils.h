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

#include "antioch_config.h"

#ifdef ANTIOCH_HAVE_EIGEN

#include <Eigen/Dense>

#include "antioch/metaprogramming.h"

namespace std {

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline
Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
max (const Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& a,
     const Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& b)
{
  using std::max;

  Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> out = a;
  const size_t size = a.size();
  for (size_t i=0; i != size; ++i)
    out[i] = max(a[i], b[i]);
  return out;
}


template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline
Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
max (const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& a,
     const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& b)
{
  using std::max;

  Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> out = a;
  const size_t size = a.size();
  for (size_t i=0; i != size; ++i)
    out[i] = max(a[i], b[i]);
  return out;
}



}


namespace Antioch
{

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline
_Scalar
max (const Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& in)
{
  return in.maxCoeff();
}

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
struct value_type<Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> >
{
  typedef Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> 
    container_type;
  typedef _Scalar type;

  static inline
  container_type
  constant(const type& in) { return container_type::Constant(in); }
};

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline
Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
zero_clone(const Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& example)
{
  if (example.size())
    return 
      Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
      (example.rows(), example.cols()).setConstant(zero_clone(example[0]));

  return 
    Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
    (example.rows(), example.cols());
}

// A function for zero-setting vectorized numeric types
template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline
void set_zero(Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& a)
{
  if (a.size())
    a = 
      value_type<Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> >::constant
        (zero_clone(a[0]));
}


template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline
_Scalar
max (const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& in)
{
  return in.maxCoeff();
}

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
struct value_type<Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> >
{
  typedef Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> 
    container_type;
  typedef _Scalar type;

  static inline
  container_type
  constant(const type& in) { return container_type::Constant(in); }
};

template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline
Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
zero_clone(const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& example)
{
  if (example.size())
    return 
      Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
      (example.rows(), example.cols()).setConstant(zero_clone(example[0]));

  return 
    Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>
    (example.rows(), example.cols());
}

// A function for zero-setting vectorized numeric types
template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline
void set_zero(Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& a)
{
  if (a.size())
    a = 
      value_type<Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> >::constant
        (zero_clone(a[0]));
}


} // end namespace Antioch

#endif // ANTIOCH_HAVE_EIGEN

#endif //ANTIOCH_EIGEN_UTILS_H
