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

#include "antioch_config.h"

#ifdef ANTIOCH_HAVE_CPPUNIT

// C++
#include <valarray>

// Antioch
#include "antioch/valarray_utils_decl.h"
#include "antioch/physical_constants.h"
#include "antioch/units.h"
#include "arrhenius_rate_vector_test_base.h"
#include "antioch/valarray_utils.h"

namespace AntiochTesting
{
  template <typename Scalar>
  class ArrheniusRateValarrayTest : public ArrheniusRateVectorTestBase<std::valarray<Scalar> >
  {
  public:

    virtual void setUp()
    {
      this->init();
      this->_example = new std::valarray<Scalar>(2*ANTIOCH_N_TUPLES);
    }

    virtual void tearDown()
    {
      this->clear();
      delete this->_example;
    }

  };

  class ArrheniusRateValarrayFloatTest : public ArrheniusRateValarrayTest<float>
  {
  public:

    CPPUNIT_TEST_SUITE( ArrheniusRateValarrayFloatTest );

    CPPUNIT_TEST( test_standard_rate );
    CPPUNIT_TEST( test_standard_deriv );
    CPPUNIT_TEST( test_standard_rate_and_deriv );

    CPPUNIT_TEST_SUITE_END();
  };

  class ArrheniusRateValarrayDoubleTest : public ArrheniusRateValarrayTest<double>
  {
  public:

    CPPUNIT_TEST_SUITE( ArrheniusRateValarrayDoubleTest );

    CPPUNIT_TEST( test_standard_rate );
    CPPUNIT_TEST( test_standard_deriv );
    CPPUNIT_TEST( test_standard_rate_and_deriv );

    CPPUNIT_TEST_SUITE_END();
  };

  class ArrheniusRateValarrayLongDoubleTest : public ArrheniusRateValarrayTest<long double>
  {
  public:

    CPPUNIT_TEST_SUITE( ArrheniusRateValarrayLongDoubleTest );

    CPPUNIT_TEST( test_standard_rate );
    CPPUNIT_TEST( test_standard_deriv );
    CPPUNIT_TEST( test_standard_rate_and_deriv );

    CPPUNIT_TEST_SUITE_END();
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( ArrheniusRateValarrayFloatTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( ArrheniusRateValarrayDoubleTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( ArrheniusRateValarrayLongDoubleTest );

} // end namespace AntiochTesting

#endif // ANTIOCH_HAVE_CPPUNIT
