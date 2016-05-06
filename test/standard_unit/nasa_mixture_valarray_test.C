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

// CppUnit
#include <cppunit/TestCase.h>

// C++
#include <valarray>

// Antioch
#include "antioch/valarray_utils_decl.h"
#include "antioch/physical_constants.h"
#include "antioch/units.h"
#include "nasa7_mixture_vector_test_base.h"
#include "nasa9_mixture_vector_test_base.h"
#include "antioch/valarray_utils.h"

namespace AntiochTesting
{
  template <typename Scalar>
  class NASA7MixtureValarrayTest : public NASA7MixtureVectorTestBase<std::valarray<Scalar> >,
                                   public CppUnit::TestCase
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

  class NASA7MixtureValarrayFloatTest : public NASA7MixtureValarrayTest<float>
  {
  public:

    CPPUNIT_TEST_SUITE( NASA7MixtureValarrayFloatTest );

    CPPUNIT_TEST( test_all_cp );
    CPPUNIT_TEST( test_all_h );
    CPPUNIT_TEST( test_all_s );

    CPPUNIT_TEST_SUITE_END();
  };

  class NASA7MixtureValarrayDoubleTest : public NASA7MixtureValarrayTest<double>
  {
  public:

    CPPUNIT_TEST_SUITE( NASA7MixtureValarrayDoubleTest );

    CPPUNIT_TEST( test_all_cp );
    CPPUNIT_TEST( test_all_h );
    CPPUNIT_TEST( test_all_s );

    CPPUNIT_TEST_SUITE_END();
  };

  class NASA7MixtureValarrayLongDoubleTest : public NASA7MixtureValarrayTest<long double>
  {
  public:

    CPPUNIT_TEST_SUITE( NASA7MixtureValarrayLongDoubleTest );

    CPPUNIT_TEST( test_all_cp );
    CPPUNIT_TEST( test_all_h );
    CPPUNIT_TEST( test_all_s );

    CPPUNIT_TEST_SUITE_END();
  };



  template <typename Scalar>
  class NASA9MixtureValarrayTest : public NASA9MixtureVectorTestBase<std::valarray<Scalar> >,
                                   public CppUnit::TestCase
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

  class NASA9MixtureValarrayFloatTest : public NASA9MixtureValarrayTest<float>
  {
  public:

    CPPUNIT_TEST_SUITE( NASA9MixtureValarrayFloatTest );

    CPPUNIT_TEST( test_all_cp );
    CPPUNIT_TEST( test_all_h );
    CPPUNIT_TEST( test_all_s );

    CPPUNIT_TEST_SUITE_END();
  };

  class NASA9MixtureValarrayDoubleTest : public NASA9MixtureValarrayTest<double>
  {
  public:

    CPPUNIT_TEST_SUITE( NASA9MixtureValarrayDoubleTest );

    CPPUNIT_TEST( test_all_cp );
    CPPUNIT_TEST( test_all_h );
    CPPUNIT_TEST( test_all_s );

    CPPUNIT_TEST_SUITE_END();
  };

  class NASA9MixtureValarrayLongDoubleTest : public NASA9MixtureValarrayTest<long double>
  {
  public:

    CPPUNIT_TEST_SUITE( NASA9MixtureValarrayLongDoubleTest );

    CPPUNIT_TEST( test_all_cp );
    CPPUNIT_TEST( test_all_h );
    CPPUNIT_TEST( test_all_s );

    CPPUNIT_TEST_SUITE_END();
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( NASA7MixtureValarrayFloatTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( NASA7MixtureValarrayDoubleTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( NASA7MixtureValarrayLongDoubleTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( NASA9MixtureValarrayFloatTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( NASA9MixtureValarrayDoubleTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( NASA9MixtureValarrayLongDoubleTest );

} // end namespace AntiochTesting

#endif // ANTIOCH_HAVE_CPPUNIT
