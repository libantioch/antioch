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

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

template <typename ReactionRate, typename Scalar>
class ReactionRateBaseTest : public CppUnit::TestCase
{
public:
  void test_rate( const ReactionRate& reaction_rate,
                  Scalar tol )
  {
    for(Scalar T = 300.1; T <= 2500.1; T += 10.)
      {
        Scalar rate = reaction_rate(T);
        Scalar exact_rate = this->exact_rate(T);

        CPPUNIT_ASSERT_DOUBLES_EQUAL( rate,
                                      exact_rate,
                                      tol );
      }
  }

  void test_deriv( const ReactionRate& reaction_rate,
                   Scalar tol )
  {
    for(Scalar T = 300.1; T <= 2500.1; T += 10.)
      {
        Scalar deriv = reaction_rate.derivative(T);
        Scalar exact_deriv = this->exact_deriv(T);

        CPPUNIT_ASSERT_DOUBLES_EQUAL( deriv,
                                      exact_deriv,
                                      tol );
      }
  }

  void test_rate_and_deriv( const ReactionRate& reaction_rate,
                            Scalar tol )
  {
    for(Scalar T = 300.1; T <= 2500.1; T += 10.)
      {
        Scalar rate;
        Scalar deriv;

        reaction_rate.rate_and_derivative(T,rate,deriv);

        Scalar exact_rate = this->exact_rate(T);
        Scalar exact_deriv = this->exact_deriv(T);

        CPPUNIT_ASSERT_DOUBLES_EQUAL( rate,
                                      exact_rate,
                                      tol );

        CPPUNIT_ASSERT_DOUBLES_EQUAL( deriv,
                                      exact_deriv,
                                      tol );
      }
  }

protected:

  virtual Scalar exact_rate( Scalar T ) =0;
  virtual Scalar exact_deriv( Scalar T ) =0;

};

#endif // ANTIOCH_HAVE_CPPUNIT
