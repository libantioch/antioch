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
#include "antioch/nasa9_curve_fit.h"
#include "antioch/temp_cache.h"
#include "nasa9_thermo_test_base.h"

namespace AntiochTesting
{
  template<typename Scalar>
  class NASA9CurveFitTest : public NASA9ThermoTestBase<Scalar>
  {
  public:

    void test_nasa9_default_temp_intervals()
    {
      {
        Antioch::NASA9CurveFit<Scalar> curve_fit( this->_all_N2_coeffs );
        this->three_interval_test(curve_fit);
      }

      {
        Antioch::NASA9CurveFit<Scalar> curve_fit( this->_all_NO2_coeffs );
        this->two_interval_test(curve_fit);
      }

    }

    void test_nasa9_user_specified_temp_intervals()
    {
      {
        std::vector<Scalar> temp(4);
        temp[0] = 200;
        temp[1] = 1000;
        temp[2] = 6000;
        temp[3] = 20000;

        Antioch::NASA9CurveFit<Scalar> curve_fit( this->_all_N2_coeffs, temp );
        this->three_interval_test(curve_fit);
      }

      {
        std::vector<Scalar> temp(3);
        temp[0] = 200;
        temp[1] = 1000;
        temp[2] = 6000;

        Antioch::NASA9CurveFit<Scalar> curve_fit( this->_all_NO2_coeffs, temp );
        this->two_interval_test(curve_fit);
      }

    }

    void three_interval_test( const Antioch::NASA9CurveFit<Scalar>& curve_fit )
    {
      Scalar T = 302.0;
      this->test_cp( T, this->_N2_coeffs_200_1000, curve_fit );
      this->test_h( T, this->_N2_coeffs_200_1000, curve_fit );
      this->test_s( T, this->_N2_coeffs_200_1000, curve_fit );

      T = 3002.0;
      this->test_cp( T, this->_N2_coeffs_1000_6000, curve_fit );
      this->test_h( T, this->_N2_coeffs_1000_6000, curve_fit );
      this->test_s( T, this->_N2_coeffs_1000_6000, curve_fit );

      T = 7002.0;
      this->test_cp( T, this->_N2_coeffs_6000_20000, curve_fit );
      this->test_h( T, this->_N2_coeffs_6000_20000, curve_fit );
      this->test_s( T, this->_N2_coeffs_6000_20000, curve_fit );
    }

    void two_interval_test( const Antioch::NASA9CurveFit<Scalar>& curve_fit )
    {
      Scalar T = 302.0;
      this->test_cp( T, this->_NO2_coeffs_200_1000, curve_fit );
      this->test_h( T, this->_NO2_coeffs_200_1000, curve_fit );
      this->test_s( T, this->_NO2_coeffs_200_1000, curve_fit );

      T = 3002.0;
      this->test_cp( T, this->_NO2_coeffs_1000_6000, curve_fit );
      this->test_h( T, this->_NO2_coeffs_1000_6000, curve_fit );
      this->test_s( T, this->_NO2_coeffs_1000_6000, curve_fit );
    }

    void test_cp( Scalar T,
                  const std::vector<Scalar>& exact_coeffs,
                 const Antioch::NASA9CurveFit<Scalar>& curve_fit )
    {
      Antioch::TempCache<Scalar> cache(T);
      Scalar cp = this->cp_exact( T,
                                  exact_coeffs[0],
                                  exact_coeffs[1],
                                  exact_coeffs[2],
                                  exact_coeffs[3],
                                  exact_coeffs[4],
                                  exact_coeffs[5],
                                  exact_coeffs[6] );

      CPPUNIT_ASSERT_DOUBLES_EQUAL( cp, curve_fit.cp_over_R(cache), this->tol() );
    }

    void test_h( Scalar T,
                 const std::vector<Scalar>& exact_coeffs,
                 const Antioch::NASA9CurveFit<Scalar>& curve_fit )
    {
      Antioch::TempCache<Scalar> cache(T);
      Scalar h = this->h_exact( T,
                                exact_coeffs[0],
                                exact_coeffs[1],
                                exact_coeffs[2],
                                exact_coeffs[3],
                                exact_coeffs[4],
                                exact_coeffs[5],
                                exact_coeffs[6],
                                exact_coeffs[7] );

      CPPUNIT_ASSERT_DOUBLES_EQUAL( h, curve_fit.h_over_RT(cache), this->tol() );
    }

    void test_s( Scalar T,
                 const std::vector<Scalar>& exact_coeffs,
                 const Antioch::NASA9CurveFit<Scalar>& curve_fit )
    {
      Antioch::TempCache<Scalar> cache(T);
      Scalar s = this->s_exact( T,
                                exact_coeffs[0],
                                exact_coeffs[1],
                                exact_coeffs[2],
                                exact_coeffs[3],
                                exact_coeffs[4],
                                exact_coeffs[5],
                                exact_coeffs[6],
                                exact_coeffs[8] );

      CPPUNIT_ASSERT_DOUBLES_EQUAL( s, curve_fit.s_over_R(cache), this->tol() );
    }

    Scalar tol()
    { return std::numeric_limits<Scalar>::epsilon() * 500; }

    Scalar cp_exact( Scalar T, Scalar a0, Scalar a1, Scalar a2,
                     Scalar a3, Scalar a4, Scalar a5, Scalar a6 )
    {
      return a0/(T*T) + a1/T + a2 + a3*T + a4*(T*T) + a5*(T*T*T) + a6*(T*T*T*T);
    }

    Scalar h_exact( Scalar T, Scalar a0, Scalar a1, Scalar a2,
                    Scalar a3, Scalar a4, Scalar a5, Scalar a6,
                    Scalar a7 )
    {
      return -a0/(T*T) + a1*std::log(T)/T + a2 + a3*T/2.0L + a4*(T*T)/3.0L
        + a5*(T*T*T)/4.0L + a6*(T*T*T*T)/5.0L + a7/T;
    }

    Scalar s_exact( Scalar T, Scalar a0, Scalar a1, Scalar a2,
                    Scalar a3, Scalar a4, Scalar a5, Scalar a6,
                    Scalar a8 )
    {
      return -a0/(2.L*T*T) - a1/T + a2*std::log(T) + a3*T + a4*(T*T)/2.0L
        + a5*(T*T*T)/3.0L + a6*(T*T*T*T)/4.0L + a8;
    }

  };

  class NASA9CurveFitTestFloat : public NASA9CurveFitTest<float>
  {
  public:
    CPPUNIT_TEST_SUITE( NASA9CurveFitTestFloat );

    CPPUNIT_TEST(test_nasa9_default_temp_intervals);
    CPPUNIT_TEST(test_nasa9_user_specified_temp_intervals);

    CPPUNIT_TEST_SUITE_END();

  };

  class NASA9CurveFitTestDouble : public NASA9CurveFitTest<double>
  {
  public:
    CPPUNIT_TEST_SUITE( NASA9CurveFitTestDouble );

    CPPUNIT_TEST(test_nasa9_default_temp_intervals);
    CPPUNIT_TEST(test_nasa9_user_specified_temp_intervals);

    CPPUNIT_TEST_SUITE_END();

  };

  class NASA9CurveFitTestLongDouble : public NASA9CurveFitTest<long double>
  {
  public:
    CPPUNIT_TEST_SUITE( NASA9CurveFitTestLongDouble );

    CPPUNIT_TEST(test_nasa9_default_temp_intervals);
    CPPUNIT_TEST(test_nasa9_user_specified_temp_intervals);

    CPPUNIT_TEST_SUITE_END();

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( NASA9CurveFitTestFloat );
  CPPUNIT_TEST_SUITE_REGISTRATION( NASA9CurveFitTestDouble );
  CPPUNIT_TEST_SUITE_REGISTRATION( NASA9CurveFitTestLongDouble );

} // end namespace AntiochTesting

#endif // ANTIOCH_HAVE_CPPUNIT
