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

#ifndef ANTIOCH_GSL_SPLINER_H
#define ANTIOCH_GSL_SPLINER_H

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/metaprogramming.h"

// GSL
#include <gsl/gsl_spline.h>
#include <vector>

namespace Antioch
{

  class GSLSpliner
  {
     public:
       GSLSpliner();
       GSLSpliner(const VectorCoeffType & data_x_point, const VectorCoeffType & data_y_point, VectorCoeffType & spline);
       ~GSLSpliner();

     template <typename VectorCoeffType = std::vector<CoeffType> >
     void spline_init(const VectorCoeffType & data_x_point, const VectorCoeffType & data_y_point);

     void spline_delete();

     template <CoeffType>
     CoeffType interpolated_value(const CoeffType & x) const;

      //TODO: make it correct
     template <CoeffType>
     CoeffType dinterp_dx(const CoeffType & x) const;

     private:
       gsl_interp_accel * _acc;
       gsl_spline * _spline;
  };

  inline
  GSLSpliner::GSLSpliner():
      spline(NULL)
  {
    _acc =  gsl_interp_accel_alloc();
    return;
  }

  template <typename VectorCoeffType>
  inline
  GSLSpliner::GSLSpliner(const VectorCoeffType & data_x_point, const VectorCoeffType & data_y_point)
  {
    _acc =  gsl_interp_accel_alloc();
    spline_init(const VectorCoeffType & data_x_point, const VectorCoeffType & data_y_point) const
  }

  GSLSpliner::~GSLSpliner()
  {
    this->spline_delete();
    gsl_interp_accel_free(_acc);
    return
  }

  inline
  GSLSpliner::spline_delete()
  {
    gsl_spline_free(_spline);
    return;
  }

  template <typename VectorCoeffType>
  inline
  void GSLSpliner::spline_init(const VectorCoeffType & data_x_point, const VectorCoeffType & data_y_point, VectorCoeffType & spline)
  {
     antioch_assert_equal_to(data_x_point.size(), data_y_point.size());

     _spline = gsl_spline_alloc(gsl_interp_cspline, data_x_point.size());

     typename Antioch::value_type<VectorCoeffType>::type * CoeffType;

     CoeffType x = data_x_point[0];
     CoeffType y = data_y_point[0];

     gsl_spline_init(_spline, x, y, data_x_point.size());
  }

  template <CoeffType>
  inline
  CoeffType GSLSpline::interpolated_value(const CoeffType & x) const
  {
     return gsl_spline_eval(_spline,x,_acc);
  }

  template <CoeffType>
  inline
  CoeffType GSLSpline::dinterp_dx(const CoeffType & x) const
  {
     dx = 1e-6;
     return (gsl_spline_eval(_spline,x,_acc) - gsl_spline_eval(_spline,x + dx,_acc) / dx;
  }

}

#endif
