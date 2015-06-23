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

#include "antioch/gsl_spliner_impl.h"

namespace Antioch
{
  namespace AntiochPrivate
  {
    GSLSplinerImplementation::GSLSplinerImplementation()
      :_acc(gsl_interp_accel_alloc()),
       _spline(NULL)
    {}

    GSLSplinerImplementation::~GSLSplinerImplementation()
    {
      this->spline_clear();
      gsl_interp_accel_free(_acc);
    }

    void GSLSplinerImplementation::spline_init( const double* x, const double* y, unsigned int size )
    {
      _spline = gsl_spline_alloc(gsl_interp_cspline, size);
      gsl_spline_init(_spline, x, y, size);
    }

    void GSLSplinerImplementation::spline_clear()
    {
      gsl_spline_free(_spline);
    }

    double GSLSplinerImplementation::eval(double x) const
    {
      return gsl_spline_eval(_spline,x,_acc);
    }

    double GSLSplinerImplementation::eval_deriv(double x) const
    {
      return gsl_spline_eval_deriv(_spline,x,_acc);
    }

  } // end namespace AntiochPrivate
} // end namespace Antioch

#endif// if ANTIOCH_HAVE_GSL
