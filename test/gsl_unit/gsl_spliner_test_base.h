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

#ifndef ANTIOCH_GSL_SPLINER_TEST_BASE_H
#define ANTIOCH_GSL_SPLINER_TEST_BASE_H

#include "antioch_config.h"

#ifdef ANTIOCH_HAVE_CPPUNIT

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

template<typename Scalar>
struct GSLSplinerTestFunction
{
  virtual Scalar operator()( const Scalar x ) =0;
};

template<typename Scalar>
struct ConstantTestFunction : public GSLSplinerTestFunction<Scalar>
{
  virtual Scalar operator()( const Scalar x )
  {
    // Scalar may actually be a PairScalar, so do this to handle those cases too.
    return Antioch::constant_clone(x,10);
  }
};

template<typename Scalar>
struct LinearTestFunction : public GSLSplinerTestFunction<Scalar>
{
  virtual Scalar operator()( const Scalar x )
  {
    // Scalar may actually be a PairScalar, so do this to handle those cases too.
    Scalar ten  = Antioch::constant_clone(x,10);
    Scalar five = Antioch::constant_clone(x,5);

    return ten + five * x;
  }
};

template<typename Scalar>
struct CubicTestFunction : public GSLSplinerTestFunction<Scalar>
{
  virtual Scalar operator()( const Scalar x )
  {
    // Scalar may actually be a PairScalar, so do this to handle those cases too.
    Scalar ten  = Antioch::constant_clone(x,10);
    Scalar two  = Antioch::constant_clone(x,2);
    Scalar five = Antioch::constant_clone(x,5);

    return ten + five * x + ten * x * x - two * x * x * x;
  }
};

template<typename Scalar>
class GSLSplinerTestBase : public CppUnit::TestCase
{
public:

  void init_data()
  {
    _n_data = 40;
    _n_test = 39;
    _x_ref.resize(_n_data);
    _y_ref.resize(_n_data);
    _x_min = -5.0L;
    _x_max = 8.0L;
  }

  void fill_ref(std::vector<Scalar>& x_ref, std::vector<Scalar>& y_ref,
                unsigned int n_data, const Scalar& x_min, const Scalar& x_max,
                GSLSplinerTestFunction<Scalar>& exact_func)
  {
    for(unsigned int i = 0; i < n_data; i++)
      {
        x_ref[i] = x_min + (Scalar)(i) * (x_max - x_min) / (Scalar)(n_data-1);
        y_ref[i] = exact_func(x_ref[i]);
      }
  }

protected:

  unsigned int _n_data;
  unsigned int _n_test;
  Scalar _x_min;
  Scalar _x_max;
  std::vector<Scalar> _x_ref;
  std::vector<Scalar> _y_ref;
};

#endif // ANTIOCH_HAVE_CPPUNIT

#endif // ANTIOCH_GSL_SPLINER_TEST_BASE_H
