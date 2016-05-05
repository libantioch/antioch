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
#include <limits>

// Antioch
#include "antioch/nasa7_curve_fit.h"
#include "antioch/temp_cache.h"
#include "nasa7_thermo_test_base.h"

namespace AntiochTesting
{
  template<typename Scalar>
  class NASA7CurveFitTest : public NASA7ThermoTestBase<Scalar>
  {
  public:

    void test_nasa7_default_temp_intervals()
    {
      Antioch::NASA7CurveFit<Scalar> curve_fit( this->_all_standard_N2_coeffs );
      this->N2_test(curve_fit);
    }

    void test_nasa7_user_specified_temp_intervals()
    {
      std::vector<Scalar> temp(3);
      temp[0] = 200;
      temp[1] = 1000;
      temp[2] = 3500;

      Antioch::NASA7CurveFit<Scalar> curve_fit( this->_all_standard_H2_coeffs, temp );
      this->H2_test(curve_fit);
    }

    void N2_test( const Antioch::NASA7CurveFit<Scalar>& curve_fit )
    {
      Scalar T = 352.0;
      this->test_cp( T, this->_N2_coeffs_300_1000, curve_fit );
      this->test_h( T, this->_N2_coeffs_300_1000, curve_fit );
      this->test_s( T, this->_N2_coeffs_300_1000, curve_fit );

      T = 3002.0;
      this->test_cp( T, this->_N2_coeffs_1000_5000, curve_fit );
      this->test_h( T, this->_N2_coeffs_1000_5000, curve_fit );
      this->test_s( T, this->_N2_coeffs_1000_5000, curve_fit );
    }

    void H2_test( const Antioch::NASA7CurveFit<Scalar>& curve_fit )
    {
      Scalar T = 352.0;
      this->test_cp( T, this->_H2_coeffs_200_1000, curve_fit );
      this->test_h( T, this->_H2_coeffs_200_1000, curve_fit );
      this->test_s( T, this->_H2_coeffs_200_1000, curve_fit );

      T = 2002.0;
      this->test_cp( T, this->_H2_coeffs_1000_3500, curve_fit );
      this->test_h( T, this->_H2_coeffs_1000_3500, curve_fit );
      this->test_s( T, this->_H2_coeffs_1000_3500, curve_fit );
    }

    void test_cp( Scalar T,
                  const std::vector<Scalar>& exact_coeffs,
                 const Antioch::NASA7CurveFit<Scalar>& curve_fit )
    {
      Antioch::TempCache<Scalar> cache(T);
      Scalar cp = this->cp_exact( T,
                                  exact_coeffs[0],
                                  exact_coeffs[1],
                                  exact_coeffs[2],
                                  exact_coeffs[3],
                                  exact_coeffs[4] );

      CPPUNIT_ASSERT_DOUBLES_EQUAL( cp, curve_fit.cp_over_R(cache), this->tol() );
    }

    void test_h( Scalar T,
                 const std::vector<Scalar>& exact_coeffs,
                 const Antioch::NASA7CurveFit<Scalar>& curve_fit )
    {
      Antioch::TempCache<Scalar> cache(T);
      Scalar h = this->h_exact( T,
                                exact_coeffs[0],
                                exact_coeffs[1],
                                exact_coeffs[2],
                                exact_coeffs[3],
                                exact_coeffs[4],
                                exact_coeffs[5] );

      CPPUNIT_ASSERT_DOUBLES_EQUAL( h, curve_fit.h_over_RT(cache), this->tol() );
    }

    void test_s( Scalar T,
                 const std::vector<Scalar>& exact_coeffs,
                 const Antioch::NASA7CurveFit<Scalar>& curve_fit )
    {
      Antioch::TempCache<Scalar> cache(T);
      Scalar s = this->s_exact( T,
                                exact_coeffs[0],
                                exact_coeffs[1],
                                exact_coeffs[2],
                                exact_coeffs[3],
                                exact_coeffs[4],
                                exact_coeffs[6] );

      CPPUNIT_ASSERT_DOUBLES_EQUAL( s, curve_fit.s_over_R(cache), this->tol() );
    }

    Scalar tol()
    { return std::numeric_limits<Scalar>::epsilon() * 500; }

  };

  class NASA7CurveFitTestFloat : public NASA7CurveFitTest<float>
  {
  public:
    CPPUNIT_TEST_SUITE( NASA7CurveFitTestFloat );

    CPPUNIT_TEST(test_nasa7_default_temp_intervals);
    CPPUNIT_TEST(test_nasa7_user_specified_temp_intervals);

    CPPUNIT_TEST_SUITE_END();

  };

  class NASA7CurveFitTestDouble : public NASA7CurveFitTest<double>
  {
  public:
    CPPUNIT_TEST_SUITE( NASA7CurveFitTestDouble );

    CPPUNIT_TEST(test_nasa7_default_temp_intervals);
    CPPUNIT_TEST(test_nasa7_user_specified_temp_intervals);

    CPPUNIT_TEST_SUITE_END();

  };

  class NASA7CurveFitTestLongDouble : public NASA7CurveFitTest<long double>
  {
  public:
    CPPUNIT_TEST_SUITE( NASA7CurveFitTestLongDouble );

    CPPUNIT_TEST(test_nasa7_default_temp_intervals);
    CPPUNIT_TEST(test_nasa7_user_specified_temp_intervals);

    CPPUNIT_TEST_SUITE_END();

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( NASA7CurveFitTestFloat );
  CPPUNIT_TEST_SUITE_REGISTRATION( NASA7CurveFitTestDouble );
  CPPUNIT_TEST_SUITE_REGISTRATION( NASA7CurveFitTestLongDouble );

} // end namespace AntiochTesting

#endif // ANTIOCH_HAVE_CPPUNIT
