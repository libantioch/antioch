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

#include "antioch_config.h" //only of we have GSL
#ifdef ANTIOCH_HAVE_GSL
#ifndef ANTIOCH_GSL_SPLINER_H
#define ANTIOCH_GSL_SPLINER_H

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/metaprogramming_decl.h"

// GSL
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <vector>

namespace Antioch
{

  template<bool B>
  struct GSLInterp
  {
      template <typename Scalar>
      Scalar interpolation(const Scalar & x, gsl_spline * spline, gsl_interp_accel * acc)
                {return gsl_spline_eval(spline,x,acc);}

      template <typename Scalar>
      Scalar dinterpolation(const Scalar & x, gsl_spline * spline, gsl_interp_accel * acc)
                {return gsl_spline_eval_deriv(spline,x,acc);}
  };

  template <>
  struct GSLInterp<true>
  {
      template <typename VectorScalar>
      VectorScalar interpolation(const VectorScalar & x, gsl_spline * spline, gsl_interp_accel * acc)
                {
                  VectorScalar out = zero_clone(x);
                  for(unsigned int i =0; i < x.size(); ++i)
                  {
                    out[i] = gsl_spline_eval(spline,x[i],acc);
                  }
                  return out;
                }

      template <typename VectorScalar>
      VectorScalar dinterpolation(const VectorScalar & x, gsl_spline * spline, gsl_interp_accel * acc)
                {
                  VectorScalar out = zero_clone(x);
                  for(unsigned int i =0; i < x.size(); ++i)
                  {
                    out[i] = gsl_spline_eval_deriv(spline,x[i],acc);
                  }
                  return out;
                }
  };

  class GSLSpliner
  {
     public:
       GSLSpliner();
       template <typename VectorCoeffType>
       GSLSpliner(const VectorCoeffType & data_x_point, const VectorCoeffType & data_y_point);
       ~GSLSpliner();

     template <typename VectorCoeffType>
     void spline_init(const VectorCoeffType & data_x_point, const VectorCoeffType & data_y_point);

     void spline_delete();

     template <typename StateType>
     StateType interpolated_value(const StateType & x) const;

     template <typename StateType>
     StateType dinterp_dx(const StateType & x) const;

     private:

       gsl_interp_accel * _acc;
       gsl_spline       * _spline;
  };

  inline
  GSLSpliner::GSLSpliner()
      :_acc(NULL),_spline(NULL)
  {
    _acc =  gsl_interp_accel_alloc();
    return;
  }

  template <typename VectorCoeffType>
  inline
  GSLSpliner::GSLSpliner(const VectorCoeffType & data_x_point, const VectorCoeffType & data_y_point)
        :_acc(NULL),_spline(NULL)
  {
    _acc =  gsl_interp_accel_alloc();
    spline_init(data_x_point, data_y_point);
  }

  GSLSpliner::~GSLSpliner()
  {
    this->spline_delete();
    gsl_interp_accel_free(_acc);
    return;
  }

  inline
  void GSLSpliner::spline_delete()
  {
    gsl_spline_free(_spline);
    return;
  }

  template <typename VectorCoeffType>
  inline
  void GSLSpliner::spline_init(const VectorCoeffType & data_x_point, const VectorCoeffType & data_y_point)
  {
     antioch_assert_equal_to(data_x_point.size(), data_y_point.size());

     _spline = gsl_spline_alloc(gsl_interp_cspline, data_x_point.size());

   // GLS takes only double, raaaaahhhh
     typedef typename rebind<VectorCoeffType,double>::type VectorGSLType;
     VectorGSLType gsl_x_point(data_x_point.size(),0);
     VectorGSLType gsl_y_point(data_y_point.size(),0);
     for(unsigned int i = 0; i < data_x_point.size(); i++)
     {
        gsl_x_point[i] = (const double)data_x_point[i];
        gsl_y_point[i] = (const double)data_y_point[i];
     }

     const double * x = &gsl_x_point[0];
     const double * y = &gsl_y_point[0];

     gsl_spline_init(_spline, x, y, data_x_point.size());
  }

  template <typename StateType>
  inline
  StateType GSLSpliner::interpolated_value(const StateType & x) const
  {
     return GSLInterp<has_size<StateType>::value>().interpolation(x, _spline, _acc);
  }

  template <typename StateType>
  inline
  StateType GSLSpliner::dinterp_dx(const StateType & x) const
  {
     return GSLInterp<has_size<StateType>::value>().dinterpolation(x, _spline, _acc);
  }

}

#endif // ANTIOCH_GSL_SPLINE

#endif// if ANTIOCH_HAVE_GSL
