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

#ifndef ANTIOCH_GSL_SPLINER_SHIM_H
#define ANTIOCH_GSL_SPLINER_SHIM_H

namespace Antioch
{
  //! Objects in this namespace are not a public API and are subject to change without notice
  namespace AntiochPrivate
  {
    // Forward declaration
    class GSLSplinerImplementation;

    //! Shim around GSLSplinerImplementation
    /*!
     * This is needed because of the inlined templated member functions GSLSpliner.
     * In other words, in order to a achieve a PIMPL design for the GSL spline
     * functionality, we need this shim (whose header will be include and functions
     * inlined in GSLSpliner so that we can forward declare GSLSplinerImplementation
     * and still hide the "raw" GSL from the users so that other packages don't need
     * to add GSL to their build system.
     */
    class GSLSplinerShim
    {
    public:

      //! Default constructor
      /*!
       * init method must be subsueqently called if using the default constructor
       */
      GSLSplinerShim();

      //! Destructor
      ~GSLSplinerShim();

      //! Initialize GSL spline objects with spline data
      void spline_init( const double* x, const double* y, unsigned int size );

      //! Clear GSL spline data. Must call spline_init again to evaluate.
      void spline_clear();

      //! Evaluate spline at point x
      double eval(double x) const;

      //! Evaluate spline derivative at point x
      double eval_deriv(double x) const;

    private:

      //! \todo This can be a unique_ptr when C++11 is mandatory
      GSLSplinerImplementation* _impl;
    };

  } // end namespace AntiochPrivate
} // end namespace Antioch

#endif // ANTIOCH_GSL_SPLINER_SHIM_H

#endif// if ANTIOCH_HAVE_GSL
