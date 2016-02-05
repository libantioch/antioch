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

// CppUnit
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

// Antioch
#include "antioch/nasa9_curve_fit.h"
#include "antioch/temp_cache.h"

namespace AntiochTesting
{
  template<typename Scalar>
  class NASA9CurveFitTest : public CppUnit::TestCase
  {
  public:

    void setUp()
    {
      this->init_N2_coeffs_200_1000();
      this->init_N2_coeffs_1000_6000();
      this->init_N2_coeffs_6000_20000();
      this->init_all_N2_coeffs();

      this->init_NO2_coeffs_200_1000();
      this->init_NO2_coeffs_1000_6000();
      this->init_all_NO2_coeffs();
    }

    void tearDown(){}

    void test_nasa9_default_temp_intervals()
    {
      {
        Antioch::NASA9CurveFit<Scalar> curve_fit( _all_N2_coeffs );
        this->three_interval_test(curve_fit);
      }

      {
        Antioch::NASA9CurveFit<Scalar> curve_fit( _all_NO2_coeffs );
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

        Antioch::NASA9CurveFit<Scalar> curve_fit( _all_N2_coeffs, temp );
        this->three_interval_test(curve_fit);
      }

      {
        std::vector<Scalar> temp(3);
        temp[0] = 200;
        temp[1] = 1000;
        temp[2] = 6000;

        Antioch::NASA9CurveFit<Scalar> curve_fit( _all_NO2_coeffs, temp );
        this->two_interval_test(curve_fit);
      }

    }

    void three_interval_test( const Antioch::NASA9CurveFit<Scalar>& curve_fit )
    {
      Scalar T = 302.0;
      this->test_cp( T, _N2_coeffs_200_1000, curve_fit );
      this->test_h( T, _N2_coeffs_200_1000, curve_fit );
      this->test_s( T, _N2_coeffs_200_1000, curve_fit );

      T = 3002.0;
      this->test_cp( T, _N2_coeffs_1000_6000, curve_fit );
      this->test_h( T, _N2_coeffs_1000_6000, curve_fit );
      this->test_s( T, _N2_coeffs_1000_6000, curve_fit );

      T = 7002.0;
      this->test_cp( T, _N2_coeffs_6000_20000, curve_fit );
      this->test_h( T, _N2_coeffs_6000_20000, curve_fit );
      this->test_s( T, _N2_coeffs_6000_20000, curve_fit );
    }

    void two_interval_test( const Antioch::NASA9CurveFit<Scalar>& curve_fit )
    {
      Scalar T = 302.0;
      this->test_cp( T, _NO2_coeffs_200_1000, curve_fit );
      this->test_h( T, _NO2_coeffs_200_1000, curve_fit );
      this->test_s( T, _NO2_coeffs_200_1000, curve_fit );

      T = 3002.0;
      this->test_cp( T, _NO2_coeffs_1000_6000, curve_fit );
      this->test_h( T, _NO2_coeffs_1000_6000, curve_fit );
      this->test_s( T, _NO2_coeffs_1000_6000, curve_fit );
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

    void init_N2_coeffs_200_1000()
    {
      _N2_coeffs_200_1000.resize(9);
      _N2_coeffs_200_1000[0] =  2.21037122e+04L;
      _N2_coeffs_200_1000[1] = -3.81846145e+02L;
      _N2_coeffs_200_1000[2] =  6.08273815e+00L;
      _N2_coeffs_200_1000[3] = -8.53091381e-03L;
      _N2_coeffs_200_1000[4] =  1.38464610e-05L;
      _N2_coeffs_200_1000[5] = -9.62579293e-09L;
      _N2_coeffs_200_1000[6] =  2.51970560e-12L;
      _N2_coeffs_200_1000[7] =  7.10845911e+02L;
      _N2_coeffs_200_1000[8] = -1.07600320e+01L;
    }

    void init_N2_coeffs_1000_6000()
    {
      _N2_coeffs_1000_6000.resize(9);
      _N2_coeffs_1000_6000[0] =  5.87709908e+05L;
      _N2_coeffs_1000_6000[1] = -2.23924255e+03L;
      _N2_coeffs_1000_6000[2] =  6.06694267e+00L;
      _N2_coeffs_1000_6000[3] = -6.13965296e-04L;
      _N2_coeffs_1000_6000[4] =  1.49179819e-07L;
      _N2_coeffs_1000_6000[5] = -1.92309442e-11L;
      _N2_coeffs_1000_6000[6] =  1.06194871e-15L;
      _N2_coeffs_1000_6000[7] =  1.28320618e+04L;
      _N2_coeffs_1000_6000[8] = -1.58663484e+01L;
    }

    void init_N2_coeffs_6000_20000()
    {
      _N2_coeffs_6000_20000.resize(9);
      _N2_coeffs_6000_20000[0] =  8.30971200e+08L;
      _N2_coeffs_6000_20000[1] = -6.42048187e+05L;
      _N2_coeffs_6000_20000[2] =  2.02020507e+02L;
      _N2_coeffs_6000_20000[3] = -3.06501961e-02L;
      _N2_coeffs_6000_20000[4] =  2.48685558e-06L;
      _N2_coeffs_6000_20000[5] = -9.70579208e-11L;
      _N2_coeffs_6000_20000[6] =  1.43751673e-15L;
      _N2_coeffs_6000_20000[7] =  4.93850663e+06L;
      _N2_coeffs_6000_20000[8] = -1.67204791e+03L;
    }

    void init_all_N2_coeffs()
    {
      _all_N2_coeffs.insert(  _all_N2_coeffs.end(),
                              _N2_coeffs_200_1000.begin(),
                              _N2_coeffs_200_1000.end() );

      _all_N2_coeffs.insert(  _all_N2_coeffs.end(),
                              _N2_coeffs_1000_6000.begin(),
                              _N2_coeffs_1000_6000.end() );

      _all_N2_coeffs.insert(  _all_N2_coeffs.end(),
                              _N2_coeffs_6000_20000.begin(),
                              _N2_coeffs_6000_20000.end() );
    }

    void init_NO2_coeffs_200_1000()
    {
      _NO2_coeffs_200_1000.resize(9);
      _NO2_coeffs_200_1000[0] = -5.64204584e+04L;
      _NO2_coeffs_200_1000[1] =  9.63309734e+02L;
      _NO2_coeffs_200_1000[2] = -2.43451851e+00L;
      _NO2_coeffs_200_1000[3] =  1.92776414e-02L;
      _NO2_coeffs_200_1000[4] = -1.87456362e-05L;
      _NO2_coeffs_200_1000[5] =  9.14553637e-09L;
      _NO2_coeffs_200_1000[6] = -1.77766146e-12L;
      _NO2_coeffs_200_1000[7] = -1.54793043e+03L;
      _NO2_coeffs_200_1000[8] =  4.06785541e+01L;
    }

    void init_NO2_coeffs_1000_6000()
    {
      _NO2_coeffs_1000_6000.resize(9);
      _NO2_coeffs_1000_6000[0] =  7.21271176e+05L;
      _NO2_coeffs_1000_6000[1] = -3.83253763e+03L;
      _NO2_coeffs_1000_6000[2] =  1.11395534e+01L;
      _NO2_coeffs_1000_6000[3] = -2.23801440e-03L;
      _NO2_coeffs_1000_6000[4] =  6.54762164e-07L;
      _NO2_coeffs_1000_6000[5] = -7.61120803e-11L;
      _NO2_coeffs_1000_6000[6] =  3.32829926e-15L;
      _NO2_coeffs_1000_6000[7] =  2.50244718e+04L;
      _NO2_coeffs_1000_6000[8] = -4.30507224e+01L;
    }

    void init_all_NO2_coeffs()
    {
      _all_NO2_coeffs.insert(  _all_NO2_coeffs.end(),
                              _NO2_coeffs_200_1000.begin(),
                              _NO2_coeffs_200_1000.end() );

      _all_NO2_coeffs.insert(  _all_NO2_coeffs.end(),
                              _NO2_coeffs_1000_6000.begin(),
                              _NO2_coeffs_1000_6000.end() );
    }

  protected:

    std::vector<Scalar> _N2_coeffs_200_1000;
    std::vector<Scalar> _N2_coeffs_1000_6000;
    std::vector<Scalar> _N2_coeffs_6000_20000;

    std::vector<Scalar> _NO2_coeffs_200_1000;
    std::vector<Scalar> _NO2_coeffs_1000_6000;

    std::vector<Scalar> _all_N2_coeffs;
    std::vector<Scalar> _all_NO2_coeffs;
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
