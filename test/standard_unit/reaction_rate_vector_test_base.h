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

#ifndef ANTIOCH_REACTION_RATE_VECTOR_TEST_BASE_H
#define ANTIOCH_REACTION_RATE_VECTOR_TEST_BASE_H

#include "antioch_config.h"

#ifdef ANTIOCH_HAVE_CPPUNIT

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

template <typename ReactionRate, typename PairScalars>
class ReactionRateVectorTestBase : public CppUnit::TestCase
{
public:

  void test_rate( const ReactionRate& reaction_rate,
                  const PairScalars& T,
                  typename Antioch::value_type<PairScalars>::type /*Scalar*/ tol )
  {
    const PairScalars rate = reaction_rate(T);
    const PairScalars exact_rate = this->exact_rate(T);

    for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
      {
        CPPUNIT_ASSERT_DOUBLES_EQUAL( rate[2*tuple],
                                      exact_rate[0],
                                      tol );

        CPPUNIT_ASSERT_DOUBLES_EQUAL( rate[2*tuple+1],
                                      exact_rate[1],
                                      tol );

      }
  }

  void test_deriv( const ReactionRate& reaction_rate,
                  const PairScalars& T,
                  typename Antioch::value_type<PairScalars>::type /*Scalar*/ tol )
  {
    const PairScalars deriv = reaction_rate.derivative(T);
    const PairScalars exact_deriv = this->exact_deriv(T);

    for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
      {
        CPPUNIT_ASSERT_DOUBLES_EQUAL( deriv[2*tuple],
                                      exact_deriv[0],
                                      tol );

        CPPUNIT_ASSERT_DOUBLES_EQUAL( deriv[2*tuple+1],
                                      exact_deriv[1],
                                      tol );

      }
  }

  void test_rate_and_deriv( const ReactionRate& reaction_rate,
                            const PairScalars& T,
                            typename Antioch::value_type<PairScalars>::type /*Scalar*/ tol )
  {
    // Init using T as the "example"
    PairScalars rate = T;
    PairScalars deriv = T;

    reaction_rate.rate_and_derivative(T,rate,deriv);

    const PairScalars exact_rate = this->exact_rate(T);
    const PairScalars exact_deriv = this->exact_deriv(T);

    for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
      {
        CPPUNIT_ASSERT_DOUBLES_EQUAL( rate[2*tuple],
                                      exact_rate[0],
                                      tol );

        CPPUNIT_ASSERT_DOUBLES_EQUAL( rate[2*tuple+1],
                                      exact_rate[1],
                                      tol );

        CPPUNIT_ASSERT_DOUBLES_EQUAL( deriv[2*tuple],
                                      exact_deriv[0],
                                      tol );

        CPPUNIT_ASSERT_DOUBLES_EQUAL( deriv[2*tuple+1],
                                      exact_deriv[1],
                                      tol );
      }
  }


protected:

  // Should be new'd/deleted in subclasses for each PairScalar type
  PairScalars* _example;

  PairScalars setup_T( const PairScalars& example )
  {
    // Construct from example to avoid resizing issues
    PairScalars T = example;
    for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
      {
        T[2*tuple]   = 1500.1;
        T[2*tuple+1] = 1600.1;
      }
    return T;
  }

  virtual PairScalars exact_rate( PairScalars T ) =0;
  virtual PairScalars exact_deriv( PairScalars T ) =0;

};

#endif // ANTIOCH_HAVE_CPPUNIT

#endif // ANTIOCH_REACTION_RATE_VECTOR_TEST_BASE_H
