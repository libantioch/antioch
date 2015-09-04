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

#ifdef ANTIOCH_HAVE_EIGEN
#include "Eigen/Dense"
#endif

#ifdef ANTIOCH_HAVE_METAPHYSICL
#include "metaphysicl/numberarray.h"
#endif

#ifdef ANTIOCH_HAVE_VEXCL
#include "vexcl/vexcl.hpp"
#endif

// Antioch
#include "antioch/eigen_utils_decl.h"
#include "antioch/metaphysicl_utils_decl.h"
#include "antioch/vexcl_utils_decl.h"
#include "antioch/valarray_utils_decl.h"
#include "antioch/vector_utils_decl.h"
#include "antioch/gsl_spliner.h"
#include "antioch/vector_utils.h"
#include "antioch/eigen_utils.h"
#include "antioch/metaphysicl_utils.h"
#include "antioch/vexcl_utils.h"
#include "antioch/valarray_utils.h"

// Base class
#include "gsl_spliner_test_base.h"

namespace AntiochTesting
{
  template<typename PairScalars>
  class GSLSplinerVecTest : public GSLSplinerTestBase<typename Antioch::value_type<PairScalars>::type>
  {
  public:

    typedef typename Antioch::value_type<PairScalars>::type Scalar;

    void init_vec_data()
    {
      this->init_data();
    }

    // Helper function
    template<typename VecFunctionType, typename ScalarFunctionType>
    void run_manually_inited_test( Scalar tol )
    {
      ScalarFunctionType exact_func;
      exact_func.init(this->_x_min, this->_x_max);

      this->fill_ref(this->_x_ref,this->_y_ref,this->_n_data, this->_x_min, this->_x_max, exact_func );

      Antioch::GSLSpliner spline;
      spline.spline_init(this->_x_ref, this->_y_ref);

      VecFunctionType vec_exact_func;
      vec_exact_func.init(this->_x_min, this->_x_max);

      this->compare_values( tol, spline, vec_exact_func );
    }

    // Helper function
    template<typename VecFunctionType, typename ScalarFunctionType>
    void run_constructor_inited_test( Scalar tol )
    {
      ScalarFunctionType exact_func;
      exact_func.init(this->_x_min, this->_x_max);

      this->fill_ref(this->_x_ref,this->_y_ref,this->_n_data, this->_x_min, this->_x_max, exact_func );

      Antioch::GSLSpliner spline(this->_x_ref, this->_y_ref);

      VecFunctionType vec_exact_func;
      vec_exact_func.init(this->_x_min, this->_x_max);

      this->compare_values( tol, spline, vec_exact_func );
    }

    // Helper function
    void compare_values( Scalar tol,
                         Antioch::GSLSpliner& spline,
                         GSLSplinerTestFunction<PairScalars>& exact_func)
    {
      // Construct from example to avoid resizing issues
      PairScalars x = *(this->_example);

      for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
        {
          x[2*tuple]   = -3.5;
          x[2*tuple+1] = 5.1;
        }

      const PairScalars gsl_value = spline.interpolated_value(x);

      const PairScalars exact_value = exact_func(x);

      for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL( gsl_value[2*tuple],
                                        exact_value[2*tuple],
                                        tol );
        }
    }

    void test_manually_inited_spline_constant_func()
    {
      const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 10;

      this->run_manually_inited_test<ConstantTestFunction<PairScalars>,ConstantTestFunction<Scalar> >(tol);
    }

    void test_constructor_inited_spline_constant_func()
    {
      const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 10;

      this->run_constructor_inited_test<ConstantTestFunction<PairScalars>,ConstantTestFunction<Scalar> >(tol);
    }

    void test_manually_inited_spline_linear_func()
    {
      const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 50;

      this->run_manually_inited_test<LinearTestFunction<PairScalars>,LinearTestFunction<Scalar> >(tol);
    }

    void test_constructor_inited_spline_linear_func()
    {
      const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 50;

      this->run_constructor_inited_test<LinearTestFunction<PairScalars>,LinearTestFunction<Scalar> >(tol);
    }

    void test_manually_inited_spline_cubic_func()
    {
      const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 50;

      this->run_manually_inited_test<CubicTestFunction<PairScalars>,CubicTestFunction<Scalar> >(tol);
    }

    void test_constructor_inited_spline_cubic_func()
    {
      const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 50;

      this->run_constructor_inited_test<CubicTestFunction<PairScalars>,CubicTestFunction<Scalar> >(tol);
    }

  protected:
    // Should be new'd/deleted in subclasses for each PairScalar type
    PairScalars* _example;
  };

  //----------------------------------------------------------------------
  // valarray tests
  //----------------------------------------------------------------------
  template <typename Scalar>
  class GslSplinerValarrayTest : public GSLSplinerVecTest<std::valarray<Scalar> >
  {
  public:

    virtual void setUp()
    {
      this->init_vec_data();
      this->_example = new std::valarray<Scalar>(2*ANTIOCH_N_TUPLES);
    }

    virtual void tearDown()
    {
      delete this->_example;
    }

  };

  class GslSplinerValarrayFloatTest : public GslSplinerValarrayTest<float>
  {
    CPPUNIT_TEST_SUITE( GslSplinerValarrayFloatTest );

    CPPUNIT_TEST( test_manually_inited_spline_constant_func );
    CPPUNIT_TEST( test_constructor_inited_spline_constant_func );
    CPPUNIT_TEST( test_manually_inited_spline_linear_func );
    CPPUNIT_TEST( test_constructor_inited_spline_linear_func );
    CPPUNIT_TEST( test_manually_inited_spline_cubic_func );
    CPPUNIT_TEST( test_constructor_inited_spline_cubic_func );

    CPPUNIT_TEST_SUITE_END();
  };

  class GslSplinerValarrayDoubleTest : public GslSplinerValarrayTest<double>
  {
    CPPUNIT_TEST_SUITE( GslSplinerValarrayDoubleTest );

    CPPUNIT_TEST( test_manually_inited_spline_constant_func );
    CPPUNIT_TEST( test_constructor_inited_spline_constant_func );
    CPPUNIT_TEST( test_manually_inited_spline_linear_func );
    CPPUNIT_TEST( test_constructor_inited_spline_linear_func );
    CPPUNIT_TEST( test_manually_inited_spline_cubic_func );
    CPPUNIT_TEST( test_constructor_inited_spline_cubic_func );

    CPPUNIT_TEST_SUITE_END();
  };

  class GslSplinerValarrayLongDoubleTest : public GslSplinerValarrayTest<long double>
  {
    CPPUNIT_TEST_SUITE( GslSplinerValarrayLongDoubleTest );

    CPPUNIT_TEST( test_manually_inited_spline_constant_func );
    CPPUNIT_TEST( test_constructor_inited_spline_constant_func );
    CPPUNIT_TEST( test_manually_inited_spline_linear_func );
    CPPUNIT_TEST( test_constructor_inited_spline_linear_func );
    CPPUNIT_TEST( test_manually_inited_spline_cubic_func );
    CPPUNIT_TEST( test_constructor_inited_spline_cubic_func );

    CPPUNIT_TEST_SUITE_END();
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( GslSplinerValarrayFloatTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( GslSplinerValarrayDoubleTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( GslSplinerValarrayLongDoubleTest );

  //----------------------------------------------------------------------
  // Eigen tests
  //----------------------------------------------------------------------
#ifdef ANTIOCH_HAVE_EIGEN

  template <typename Scalar>
  class GslSplinerEigenTest : public GSLSplinerVecTest<Eigen::Array<Scalar,2*ANTIOCH_N_TUPLES,1> >
  {
  public:
    virtual void setUp()
    {
      this->init_vec_data();
      this->_example = new Eigen::Array<Scalar, 2*ANTIOCH_N_TUPLES, 1>();
    }

    virtual void tearDown()
    {
      delete this->_example;
    }
  };

  class GslSplinerEigenFloatTest : public GslSplinerEigenTest<float>
  {
    CPPUNIT_TEST_SUITE( GslSplinerEigenFloatTest );

    CPPUNIT_TEST( test_manually_inited_spline_constant_func );
    CPPUNIT_TEST( test_constructor_inited_spline_constant_func );
    CPPUNIT_TEST( test_manually_inited_spline_linear_func );
    CPPUNIT_TEST( test_constructor_inited_spline_linear_func );
    CPPUNIT_TEST( test_manually_inited_spline_cubic_func );
    CPPUNIT_TEST( test_constructor_inited_spline_cubic_func );

    CPPUNIT_TEST_SUITE_END();
  };

  class GslSplinerEigenDoubleTest : public GslSplinerEigenTest<double>
  {
    CPPUNIT_TEST_SUITE( GslSplinerEigenDoubleTest );

    CPPUNIT_TEST( test_manually_inited_spline_constant_func );
    CPPUNIT_TEST( test_constructor_inited_spline_constant_func );
    CPPUNIT_TEST( test_manually_inited_spline_linear_func );
    CPPUNIT_TEST( test_constructor_inited_spline_linear_func );
    CPPUNIT_TEST( test_manually_inited_spline_cubic_func );
    CPPUNIT_TEST( test_constructor_inited_spline_cubic_func );

    CPPUNIT_TEST_SUITE_END();
  };

  class GslSplinerEigenLongDoubleTest : public GslSplinerEigenTest<long double>
  {
    CPPUNIT_TEST_SUITE( GslSplinerEigenLongDoubleTest );

    CPPUNIT_TEST( test_manually_inited_spline_constant_func );
    CPPUNIT_TEST( test_constructor_inited_spline_constant_func );
    CPPUNIT_TEST( test_manually_inited_spline_linear_func );
    CPPUNIT_TEST( test_constructor_inited_spline_linear_func );
    CPPUNIT_TEST( test_manually_inited_spline_cubic_func );
    CPPUNIT_TEST( test_constructor_inited_spline_cubic_func );

    CPPUNIT_TEST_SUITE_END();
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( GslSplinerEigenFloatTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( GslSplinerEigenDoubleTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( GslSplinerEigenLongDoubleTest );

#endif // ANTIOCH_HAVE_EIGEN


  //----------------------------------------------------------------------
  // MetaPhysicL tests
  //----------------------------------------------------------------------
#ifdef ANTIOCH_HAVE_METAPHYSICL

  template <typename Scalar>
  class GslSplinerMetaPhysicLTest : public GSLSplinerVecTest<MetaPhysicL::NumberArray<2*ANTIOCH_N_TUPLES,Scalar> >
  {
  public:
    virtual void setUp()
    {
      this->init_vec_data();
      this->_example = new MetaPhysicL::NumberArray<2*ANTIOCH_N_TUPLES,Scalar>(0);
    }

    virtual void tearDown()
    {
      delete this->_example;
    }
  };

  class GslSplinerMetaPhysicLFloatTest : public GslSplinerMetaPhysicLTest<float>
  {
    CPPUNIT_TEST_SUITE( GslSplinerMetaPhysicLFloatTest );

    CPPUNIT_TEST( test_manually_inited_spline_constant_func );
    CPPUNIT_TEST( test_constructor_inited_spline_constant_func );
    CPPUNIT_TEST( test_manually_inited_spline_linear_func );
    CPPUNIT_TEST( test_constructor_inited_spline_linear_func );
    CPPUNIT_TEST( test_manually_inited_spline_cubic_func );
    CPPUNIT_TEST( test_constructor_inited_spline_cubic_func );

    CPPUNIT_TEST_SUITE_END();
  };

  class GslSplinerMetaPhysicLDoubleTest : public GslSplinerMetaPhysicLTest<double>
  {
    CPPUNIT_TEST_SUITE( GslSplinerMetaPhysicLDoubleTest );

    CPPUNIT_TEST( test_manually_inited_spline_constant_func );
    CPPUNIT_TEST( test_constructor_inited_spline_constant_func );
    CPPUNIT_TEST( test_manually_inited_spline_linear_func );
    CPPUNIT_TEST( test_constructor_inited_spline_linear_func );
    CPPUNIT_TEST( test_manually_inited_spline_cubic_func );
    CPPUNIT_TEST( test_constructor_inited_spline_cubic_func );

    CPPUNIT_TEST_SUITE_END();
  };

  class GslSplinerMetaPhysicLLongDoubleTest : public GslSplinerMetaPhysicLTest<long double>
  {
    CPPUNIT_TEST_SUITE( GslSplinerMetaPhysicLLongDoubleTest );

    CPPUNIT_TEST( test_manually_inited_spline_constant_func );
    CPPUNIT_TEST( test_constructor_inited_spline_constant_func );
    CPPUNIT_TEST( test_manually_inited_spline_linear_func );
    CPPUNIT_TEST( test_constructor_inited_spline_linear_func );
    CPPUNIT_TEST( test_manually_inited_spline_cubic_func );
    CPPUNIT_TEST( test_constructor_inited_spline_cubic_func );

    CPPUNIT_TEST_SUITE_END();
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( GslSplinerMetaPhysicLFloatTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( GslSplinerMetaPhysicLDoubleTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( GslSplinerMetaPhysicLLongDoubleTest );

#endif // ANTIOCH_HAVE_METAPHYSICL


  //----------------------------------------------------------------------
  // VexCL tests
  //----------------------------------------------------------------------
#ifdef ANTIOCH_HAVE_VEXCL

  class GslSplinerVexCLFloatTest : public GSLSplinerVecTest<vex::vector<float> >
  {
    CPPUNIT_TEST_SUITE( GslSplinerVexCLFloatTest );

    CPPUNIT_TEST( test_manually_inited_spline_constant_func );
    CPPUNIT_TEST( test_constructor_inited_spline_constant_func );
    CPPUNIT_TEST( test_manually_inited_spline_linear_func );
    CPPUNIT_TEST( test_constructor_inited_spline_linear_func );
    CPPUNIT_TEST( test_manually_inited_spline_cubic_func );
    CPPUNIT_TEST( test_constructor_inited_spline_cubic_func );

    CPPUNIT_TEST_SUITE_END();

  private:

    vex::Context* _ctx;

  public:

    virtual void setUp()
    {
      this->init_vec_data();
      _ctx = new vex::Context(vex::Filter::All);
      this->_example = new vex::vector<float>(*_ctx, 2*ANTIOCH_N_TUPLES);
    }

    virtual void tearDown()
    {
      delete _ctx;
      delete this->_example;
    }
  };

  class GslSplinerVexCLDoubleTest : public GSLSplinerVecTest<vex::vector<double> >
  {
    CPPUNIT_TEST_SUITE( GslSplinerVexCLDoubleTest );

    CPPUNIT_TEST( test_manually_inited_spline_constant_func );
    CPPUNIT_TEST( test_constructor_inited_spline_constant_func );
    CPPUNIT_TEST( test_manually_inited_spline_linear_func );
    CPPUNIT_TEST( test_constructor_inited_spline_linear_func );
    CPPUNIT_TEST( test_manually_inited_spline_cubic_func );
    CPPUNIT_TEST( test_constructor_inited_spline_cubic_func );

    CPPUNIT_TEST_SUITE_END();

  private:

    vex::Context* _ctx;

  public:

    virtual void setUp()
    {
      this->init_vec_data();
      _ctx = new vex::Context(vex::Filter::DoublePrecision);
      this->_example = new vex::vector<double>(*_ctx, 2*ANTIOCH_N_TUPLES);
    }

    virtual void tearDown()
    {
      delete _ctx;
      delete this->_example;
    }
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( GslSplinerVexCLFloatTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( GslSplinerVexCLDoubleTest );

#endif // ANTIOCH_HAVE_VEXCL

} // end namespace AntiochTesting

#endif // ANTIOCH_HAVE_CPPUNIT
