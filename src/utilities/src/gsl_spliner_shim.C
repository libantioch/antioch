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

#include "antioch_config.h" //only of we have GSL

#ifdef ANTIOCH_HAVE_GSL

// This class
#include "antioch/gsl_spliner_shim.h"

// Antioch
#include "antioch/gsl_spliner_impl.h"

namespace Antioch
{
  //! Objects in this namespace are not a public API and are subject to change without notice
  namespace AntiochPrivate
  {
    GSLSplinerShim::GSLSplinerShim()
      : _impl( new GSLSplinerImplementation() )
    {}

    GSLSplinerShim::~GSLSplinerShim()
    {
      delete _impl;
    }

    void GSLSplinerShim::spline_init( const double* x, const double* y, unsigned int size )
    {
      _impl->spline_init(x,y,size);
    }

    void GSLSplinerShim::spline_clear()
    {
      _impl->spline_clear();
    }

    double GSLSplinerShim::eval(double x) const
    {
      return _impl->eval(x);
    }

    double GSLSplinerShim::eval_deriv(double x) const
    {
      return _impl->eval_deriv(x);
    }

  } // end namespace AntiochPrivate
} // end namespace Antioch

#endif// if ANTIOCH_HAVE_GSL
