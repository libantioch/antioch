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

#include "antioch_config.h"

#ifdef ANTIOCH_HAVE_CPPUNIT

// C++
#include <limits>
#include <cmath>

// Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/gsl_spliner.h"
#include "antioch/vector_utils.h"

// Base class
#include "gsl_spliner_test_base.h"

namespace AntiochTesting
{
  template<typename Scalar>
  class GSLSplinerTest : public GSLSplinerTestBase<Scalar>
  {
  public:

    void setUp()
    {
      this->init_data();
    }

    // Helper function
    template<typename FunctionType>
    void run_manually_inited_test( Scalar tol )
    {
      FunctionType exact_func;
      exact_func.init(this->_x_min, this->_x_max);

      this->fill_ref(this->_x_ref,this->_y_ref,this->_n_data, this->_x_min, this->_x_max, exact_func );

      Antioch::GSLSpliner spline;
      spline.spline_init(this->_x_ref, this->_y_ref);

      this->compare_values( tol, spline, exact_func );
    }

    // Helper function
    template<typename FunctionType>
    void run_constructor_inited_test( Scalar tol )
    {
      FunctionType exact_func;
      exact_func.init(this->_x_min, this->_x_max);

      this->fill_ref(this->_x_ref,this->_y_ref,this->_n_data, this->_x_min, this->_x_max, exact_func );

      Antioch::GSLSpliner spline(this->_x_ref, this->_y_ref);

      this->compare_values( tol, spline, exact_func );
    }

    // Helper function
    void compare_values( Scalar tol, Antioch::GSLSpliner& spline, GSLSplinerTestFunction<Scalar>& exact_func )
    {
      for(unsigned int n = 0; n < this->_n_test; n++)
        {
          Scalar x = this->_x_min + (Scalar)(n) * (this->_x_max - this->_x_min) / (Scalar)(this->_n_test-1);
          Scalar exact_value = exact_func(x);
          Scalar interp_value = spline.interpolated_value(x);
          CPPUNIT_ASSERT_DOUBLES_EQUAL( interp_value,
                                        exact_value,
                                        tol );
        }
    }

    void test_manually_inited_spline_constant_func()
    {
      const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 10;

      this->run_manually_inited_test<ConstantTestFunction<Scalar> >(tol);
    }

    void test_constructor_inited_spline_constant_func()
    {
      const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 10;

      this->run_constructor_inited_test<ConstantTestFunction<Scalar> >(tol);
    }

    void test_manually_inited_spline_linear_func()
    {
      const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 50;

      this->run_manually_inited_test<LinearTestFunction<Scalar> >(tol);
    }

    void test_constructor_inited_spline_linear_func()
    {
      const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 50;

      this->run_constructor_inited_test<LinearTestFunction<Scalar> >(tol);
    }

    void test_manually_inited_spline_cubic_func()
    {
      const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 50;

      this->run_manually_inited_test<CubicTestFunction<Scalar> >(tol);
    }

    void test_constructor_inited_spline_cubic_func()
    {
      const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 50;

      this->run_constructor_inited_test<CubicTestFunction<Scalar> >(tol);
    }
  };

  class GslSplinerFloatTest : public GSLSplinerTest<float>
  {
  public:
    CPPUNIT_TEST_SUITE( GslSplinerFloatTest );

    CPPUNIT_TEST( test_manually_inited_spline_constant_func );
    CPPUNIT_TEST( test_constructor_inited_spline_constant_func );
    CPPUNIT_TEST( test_manually_inited_spline_linear_func );
    CPPUNIT_TEST( test_constructor_inited_spline_linear_func );
    CPPUNIT_TEST( test_manually_inited_spline_cubic_func );
    CPPUNIT_TEST( test_constructor_inited_spline_cubic_func );

    CPPUNIT_TEST_SUITE_END();
  };

  class GslSplinerDoubleTest : public GSLSplinerTest<double>
  {
  public:
    CPPUNIT_TEST_SUITE( GslSplinerDoubleTest );

    CPPUNIT_TEST( test_manually_inited_spline_constant_func );
    CPPUNIT_TEST( test_constructor_inited_spline_constant_func );
    CPPUNIT_TEST( test_manually_inited_spline_linear_func );
    CPPUNIT_TEST( test_constructor_inited_spline_linear_func );
    CPPUNIT_TEST( test_manually_inited_spline_cubic_func );
    CPPUNIT_TEST( test_constructor_inited_spline_cubic_func );

    CPPUNIT_TEST_SUITE_END();
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( GslSplinerFloatTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( GslSplinerDoubleTest );

} // end namespace AntiochTesting

#endif // ANTIOCH_HAVE_CPPUNIT
