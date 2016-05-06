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

#ifdef ANTIOCH_HAVE_EIGEN

// Eigen
#include "Eigen/Dense"

// CppUnit
#include <cppunit/TestCase.h>

// Antioch
#include "antioch/eigen_utils_decl.h"
#include "antioch/physical_constants.h"
#include "antioch/units.h"
#include "nasa7_mixture_vector_test_base.h"
#include "nasa9_mixture_vector_test_base.h"
#include "antioch/eigen_utils.h"

namespace AntiochTesting
{
  template <typename Scalar>
  class NASA7MixtureEigenTest : public NASA7MixtureVectorTestBase<Eigen::Array<Scalar,2*ANTIOCH_N_TUPLES,1> >,
                                   public CppUnit::TestCase
  {
  public:

    virtual void setUp()
    {
      this->init();
      this->_example = new Eigen::Array<Scalar, 2*ANTIOCH_N_TUPLES, 1>();
    }

    virtual void tearDown()
    {
      this->clear();
      delete this->_example;
    }

  };

  class NASA7MixtureEigenFloatTest : public NASA7MixtureEigenTest<float>
  {
  public:

    CPPUNIT_TEST_SUITE( NASA7MixtureEigenFloatTest );

    CPPUNIT_TEST( test_all_cp );
    CPPUNIT_TEST( test_all_h );
    CPPUNIT_TEST( test_all_s );

    CPPUNIT_TEST_SUITE_END();
  };

  class NASA7MixtureEigenDoubleTest : public NASA7MixtureEigenTest<double>
  {
  public:

    CPPUNIT_TEST_SUITE( NASA7MixtureEigenDoubleTest );

    CPPUNIT_TEST( test_all_cp );
    CPPUNIT_TEST( test_all_h );
    CPPUNIT_TEST( test_all_s );

    CPPUNIT_TEST_SUITE_END();
  };

  class NASA7MixtureEigenLongDoubleTest : public NASA7MixtureEigenTest<long double>
  {
  public:

    CPPUNIT_TEST_SUITE( NASA7MixtureEigenLongDoubleTest );

    CPPUNIT_TEST( test_all_cp );
    CPPUNIT_TEST( test_all_h );
    CPPUNIT_TEST( test_all_s );

    CPPUNIT_TEST_SUITE_END();
  };



  template <typename Scalar>
  class NASA9MixtureEigenTest : public NASA9MixtureVectorTestBase<Eigen::Array<Scalar,2*ANTIOCH_N_TUPLES,1> >,
                                   public CppUnit::TestCase
  {
  public:

    virtual void setUp()
    {
      this->init();
      this->_example = new Eigen::Array<Scalar, 2*ANTIOCH_N_TUPLES, 1>();
    }

    virtual void tearDown()
    {
      this->clear();
      delete this->_example;
    }

  };

  class NASA9MixtureEigenFloatTest : public NASA9MixtureEigenTest<float>
  {
  public:

    CPPUNIT_TEST_SUITE( NASA9MixtureEigenFloatTest );

    CPPUNIT_TEST( test_all_cp );
    CPPUNIT_TEST( test_all_h );
    CPPUNIT_TEST( test_all_s );

    CPPUNIT_TEST_SUITE_END();
  };

  class NASA9MixtureEigenDoubleTest : public NASA9MixtureEigenTest<double>
  {
  public:

    CPPUNIT_TEST_SUITE( NASA9MixtureEigenDoubleTest );

    CPPUNIT_TEST( test_all_cp );
    CPPUNIT_TEST( test_all_h );
    CPPUNIT_TEST( test_all_s );

    CPPUNIT_TEST_SUITE_END();
  };

  class NASA9MixtureEigenLongDoubleTest : public NASA9MixtureEigenTest<long double>
  {
  public:

    CPPUNIT_TEST_SUITE( NASA9MixtureEigenLongDoubleTest );

    CPPUNIT_TEST( test_all_cp );
    CPPUNIT_TEST( test_all_h );
    CPPUNIT_TEST( test_all_s );

    CPPUNIT_TEST_SUITE_END();
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( NASA7MixtureEigenFloatTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( NASA7MixtureEigenDoubleTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( NASA7MixtureEigenLongDoubleTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( NASA9MixtureEigenFloatTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( NASA9MixtureEigenDoubleTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( NASA9MixtureEigenLongDoubleTest );

} // end namespace AntiochTesting

#endif // ANTIOCH_HAVE_EIGEN

#endif // ANTIOCH_HAVE_CPPUNIT
