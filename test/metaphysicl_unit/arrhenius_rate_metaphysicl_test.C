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

#ifdef ANTIOCH_HAVE_METAPHYSICL

// MetaPhysicL
#include "metaphysicl/numberarray.h"

// Antioch
#include "antioch/metaphysicl_utils_decl.h"
#include "antioch/physical_constants.h"
#include "antioch/units.h"
#include "arrhenius_rate_vector_test_base.h"
#include "antioch/metaphysicl_utils.h"

namespace AntiochTesting
{
  template <typename Scalar>
  class ArrheniusRateMetaPhysicLTest : public ArrheniusRateVectorTestBase<MetaPhysicL::NumberArray<2*ANTIOCH_N_TUPLES, float> >
  {
  public:

    virtual void setUp()
    {
      this->init();
      this->_example = new MetaPhysicL::NumberArray<2*ANTIOCH_N_TUPLES, float>(0);
    }

    virtual void tearDown()
    {
      this->clear();
      delete this->_example;
    }

  };

  class ArrheniusRateMetaPhysicLFloatTest : public ArrheniusRateMetaPhysicLTest<float>
  {
  public:

    CPPUNIT_TEST_SUITE( ArrheniusRateMetaPhysicLFloatTest );

    CPPUNIT_TEST( test_standard_rate );
    CPPUNIT_TEST( test_standard_deriv );
    CPPUNIT_TEST( test_standard_rate_and_deriv );

    CPPUNIT_TEST_SUITE_END();
  };

  class ArrheniusRateMetaPhysicLDoubleTest : public ArrheniusRateMetaPhysicLTest<double>
  {
  public:

    CPPUNIT_TEST_SUITE( ArrheniusRateMetaPhysicLDoubleTest );

    CPPUNIT_TEST( test_standard_rate );
    CPPUNIT_TEST( test_standard_deriv );
    CPPUNIT_TEST( test_standard_rate_and_deriv );

    CPPUNIT_TEST_SUITE_END();
  };

  class ArrheniusRateMetaPhysicLLongDoubleTest : public ArrheniusRateMetaPhysicLTest<long double>
  {
  public:

    CPPUNIT_TEST_SUITE( ArrheniusRateMetaPhysicLLongDoubleTest );

    CPPUNIT_TEST( test_standard_rate );
    CPPUNIT_TEST( test_standard_deriv );
    CPPUNIT_TEST( test_standard_rate_and_deriv );

    CPPUNIT_TEST_SUITE_END();
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( ArrheniusRateMetaPhysicLFloatTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( ArrheniusRateMetaPhysicLDoubleTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( ArrheniusRateMetaPhysicLLongDoubleTest );

} // end namespace AntiochTesting

#endif // ANTIOCH_HAVE_METAPHYSICL

#endif // ANTIOCH_HAVE_CPPUNIT
