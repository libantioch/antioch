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

#ifndef ANTIOCH_GSL_SPLINER_IMPL_H
#define ANTIOCH_GSL_SPLINER_IMPL_H

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/metaprogramming_decl.h"

// GSL
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

namespace Antioch
{
  //! Objects in this namespace are not a public API and are subject to change without notice
  namespace AntiochPrivate
  {
    //! PIMPL wrapper for raw GSL objects
    class GSLSplinerImplementation
    {
    public:

      //! Default constructor
      /*!
       * init method must be subsueqently called if using the default constructor
       */
      GSLSplinerImplementation();

      //! Destructor
      ~GSLSplinerImplementation();

      //! Initialize GSL spline objects with spline data
      void spline_init( const double* x, const double* y, unsigned int size );

      //! Clear GSL spline data. Must call spline_init again to evaluate.
      void spline_clear();

      //! Evaluate spline at point x
      double eval(double x) const;

      //! Evaluate spline derivative at point x
      double eval_deriv(double x) const;

    private:

      gsl_interp_accel * _acc;
      gsl_spline       * _spline;

    };

  } // end namespace AntiochPrivate

} // end namespace Antioch

#endif // ANTIOCH_GSL_SPLINER_IMPL_H

#endif// if ANTIOCH_HAVE_GSL
