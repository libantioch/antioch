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

#ifndef ANTIOCH_ARRHENIUS_RATE_VECTOR_TEST_BASE_H
#define ANTIOCH_ARRHENIUS_RATE_VECTOR_TEST_BASE_H

#include "antioch_config.h"

#ifdef ANTIOCH_HAVE_CPPUNIT

// C++
#include <limits>

// Antioch
#include "antioch/physical_constants.h"
#include "antioch/units.h"

// Base class
#include "arrhenius_rate_test_helper.h"
#include "reaction_rate_vector_test_base.h"

template<typename PairScalars>
class ArrheniusRateVectorTestBase : public ArrheniusRateTestHelper</*Scalar*/typename Antioch::value_type<PairScalars>::type>,
                                    public ReactionRateVectorBaseTest<Antioch::ArrheniusRate</*Scalar*/typename Antioch::value_type<PairScalars>::type>,PairScalars>
{
public:
  virtual void init()
  {
    typedef typename Antioch::value_type<PairScalars>::type Scalar;

    Scalar Cf = 1.4L;
    Scalar Ea = 5.0L;
    Scalar R = 1.0L; // Ea in K

    this->reset_params(Cf,Ea,R);
    _rate = new Antioch::ArrheniusRate<Scalar>(Cf,Ea,R);
  }

  virtual void clear()
  {
    delete _rate;
  }

  void test_standard_rate()
  {
    typedef typename Antioch::value_type<PairScalars>::type Scalar;

    const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;
    PairScalars T = this->setup_T(*(this->_example));
    this->test_rate( *_rate, T, tol );
  }

  void test_standard_deriv()
  {
    typedef typename Antioch::value_type<PairScalars>::type Scalar;

    const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;
    PairScalars T = this->setup_T(*(this->_example));
    this->test_deriv( *_rate, T, tol );
  }

  void test_standard_rate_and_deriv()
  {
    typedef typename Antioch::value_type<PairScalars>::type Scalar;

    const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;
    PairScalars T = this->setup_T(*(this->_example));
    this->test_rate_and_deriv( *_rate, T, tol );
  }

protected:

  Antioch::ArrheniusRate<typename Antioch::value_type<PairScalars>::type>* _rate;

  virtual PairScalars exact_rate( PairScalars T )
  {
    PairScalars e_rate = *(this->_example);
    for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
      {
        e_rate[2*tuple]   =  this->value(T[2*tuple]);
        e_rate[2*tuple+1] =  this->value(T[2*tuple+1]);
      }
    return e_rate;
  }

  virtual PairScalars exact_deriv( PairScalars T )
  {
    PairScalars e_deriv = *(this->_example);
    for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
      {
        e_deriv[2*tuple]   =  this->deriv(T[2*tuple]);
        e_deriv[2*tuple+1] =  this->deriv(T[2*tuple+1]);
      }
    return e_deriv;
  }

};

#endif // ANTIOCH_HAVE_CPPUNIT

#endif // ANTIOCH_ARRHENIUS_RATE_VECTOR_TEST_BASE_H
