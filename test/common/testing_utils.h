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

#ifndef ANTIOCH_TESTING_UTILS_H
#define ANTIOCH_TESTING_UTILS_H

#include "antioch_config.h"

#ifdef ANTIOCH_HAVE_CPPUNIT

// CppUnit
#include <cppunit/extensions/HelperMacros.h>

// C++
#include <cmath>

namespace AntiochTesting
{
  template<typename Scalar>
  class TestingUtilities
  {
  public:

    void test_scalar_rel( Scalar expected, Scalar actual, Scalar tol )
    {
      Scalar rel_error = std::abs( (actual-expected) );

      // Don't normalize if we're expecting the result to be really small
      if( expected > tol )
        rel_error = std::abs( rel_error/expected);

      CPPUNIT_ASSERT_DOUBLES_EQUAL(Scalar(0), rel_error, tol );
    }
  };
} // end namespace AntiochTesting

#endif // ANTIOCH_HAVE_CPPUNIT

#endif // ANTIOCH_TESTING_UTILS_H
