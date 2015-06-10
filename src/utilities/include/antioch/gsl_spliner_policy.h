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

#ifndef ANTIOCH_GSL_SPLINER_POLICY_H
#define ANTIOCH_GSL_SPLINER_POLICY_H

// GSL
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

namespace Antioch
{
  namespace AntiochPrivate
  {

    template<bool B>
    struct GSLSplinerPolicy
    {
      template <typename Scalar>
      Scalar interpolation(const Scalar & x, gsl_spline * spline, gsl_interp_accel * acc)
      {return gsl_spline_eval(spline,x,acc);}

      template <typename Scalar>
      Scalar dinterpolation(const Scalar & x, gsl_spline * spline, gsl_interp_accel * acc)
      {return gsl_spline_eval_deriv(spline,x,acc);}
    };

    template <>
    struct GSLSplinerPolicy<true>
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

  } // end namespace AntiochPrivate

} // end namespace Antioch

#endif // ANTIOCH_GSL_SPLINER_POLICY_H


#endif// if ANTIOCH_HAVE_GSL
