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

// Antioch
#include "antioch/physical_constants.h"
#include "antioch/units.h"

// Base class
#include "arrhenius_rate_test_helper.h"
#include "reaction_rate_test_base.h"

using namespace Antioch;

template<typename Scalar>
class ArrheniusRateTest : public ArrheniusRateTestHelper<Scalar>,
                          public ReactionRateTestBase<Antioch::ArrheniusRate<Scalar>,Scalar>
{
public:
  void setUp()
  {
    Scalar Cf = 1.4L;
    Scalar Ea = 298.0L;
    Scalar R = 1.0L; // Ea in K

    this->reset_params(Cf,Ea,R);
    _rate = new ArrheniusRate<Scalar>(Cf,Ea,R);
  }

  void tearDown()
  {
    delete _rate;
   }

  void test_standard()
  {
    const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;

    this->test_rate( *_rate, tol );
    this->test_deriv( *_rate, tol );
    this->test_rate_and_deriv( *_rate, tol );
  }

  void test_reset_scalar_params()
  {
    Scalar Cf = 1e-7L;
    Scalar Ea = 36000.L;
    Scalar R = Constants::R_universal<Scalar>()*Units<Scalar>("cal").get_SI_factor();

    this->reset_params( Cf, Ea, R );

    _rate->set_Cf(Cf);
    _rate->set_Ea(Ea);
    _rate->set_rscale(R);

    const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;

    this->test_rate( *_rate, tol );
    this->test_deriv( *_rate, tol );
    this->test_rate_and_deriv( *_rate, tol );
  }

  void test_reset_vector_params2()
  {
    Scalar Cf = 2.5e-7L;
    Scalar Ea = 43000.L; //still in cal

    this->reset_params( Cf, Ea );

    std::vector<Scalar> values(2);
    values[0] = Cf;
    values[1] = Ea;
    _rate->reset_coefs(values);

    const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;

    this->test_rate( *_rate, tol );
    this->test_deriv( *_rate, tol );
    this->test_rate_and_deriv( *_rate, tol );
  }

  void test_reset_vector_params3()
  {
    Scalar Cf = 2.1e-11L;
    Scalar Ea = 100000.L;
    Scalar R = Constants::R_universal<Scalar>();

    this->reset_params( Cf, Ea, R );

    std::vector<Scalar> values(3);
    values[0] = Cf;
    values[1] = Ea;
    values[2] = R;
    _rate->reset_coefs(values);

    const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;

    this->test_rate( *_rate, tol );
    this->test_deriv( *_rate, tol );
    this->test_rate_and_deriv( *_rate, tol );
  }

protected:

  ArrheniusRate<Scalar>* _rate;

  virtual Scalar exact_rate( Scalar T )
  {
    return this->value(T);
  }

  virtual Scalar exact_deriv( Scalar T )
  {
    return this->deriv(T);
  }

};

class ArrheniusRateFloatTest : public ArrheniusRateTest<float>
{
public:
  CPPUNIT_TEST_SUITE( ArrheniusRateFloatTest );

  CPPUNIT_TEST( test_standard );
  CPPUNIT_TEST( test_reset_scalar_params );
  CPPUNIT_TEST( test_reset_vector_params2 );
  CPPUNIT_TEST( test_reset_vector_params3 );

  CPPUNIT_TEST_SUITE_END();
};

class ArrheniusRateDoubleTest : public ArrheniusRateTest<double>
{
public:
  CPPUNIT_TEST_SUITE( ArrheniusRateDoubleTest );

  CPPUNIT_TEST( test_standard );
  CPPUNIT_TEST( test_reset_scalar_params );
  CPPUNIT_TEST( test_reset_vector_params2 );
  CPPUNIT_TEST( test_reset_vector_params3 );

  CPPUNIT_TEST_SUITE_END();
};

class ArrheniusRateLongDoubleTest : public ArrheniusRateTest<long double>
{
public:
  CPPUNIT_TEST_SUITE( ArrheniusRateLongDoubleTest );

  CPPUNIT_TEST( test_standard );
  CPPUNIT_TEST( test_reset_scalar_params );
  CPPUNIT_TEST( test_reset_vector_params2 );
  CPPUNIT_TEST( test_reset_vector_params3 );

  CPPUNIT_TEST_SUITE_END();
};

CPPUNIT_TEST_SUITE_REGISTRATION( ArrheniusRateFloatTest );
CPPUNIT_TEST_SUITE_REGISTRATION( ArrheniusRateDoubleTest );
CPPUNIT_TEST_SUITE_REGISTRATION( ArrheniusRateLongDoubleTest );

#endif // ANTIOCH_HAVE_CPPUNIT
