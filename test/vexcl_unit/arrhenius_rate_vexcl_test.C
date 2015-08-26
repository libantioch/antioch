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

#ifdef ANTIOCH_HAVE_VEXCL

// VexCL
#include "vexcl/vexcl.hpp"

// Antioch
#include "antioch/vexcl_utils_decl.h"
#include "antioch/physical_constants.h"
#include "antioch/units.h"
#include "arrhenius_rate_vector_test_base.h"
#include "antioch/vexcl_utils.h"

class ArrheniusRateVexCLFloatTest : public ArrheniusRateVectorTestBase<vex::vector<float> >
{
public:

  CPPUNIT_TEST_SUITE( ArrheniusRateVexCLFloatTest );

  CPPUNIT_TEST( test_standard_rate );
  CPPUNIT_TEST( test_standard_deriv );
  CPPUNIT_TEST( test_standard_rate_and_deriv );

  CPPUNIT_TEST_SUITE_END();

private:

  vex::Context* _ctx;

public:

  virtual void setUp()
  {
    this->init();
    _ctx = new vex::Context(vex::Filter::All);
    this->_example = new vex::vector<float>(*_ctx, 2*ANTIOCH_N_TUPLES);
  }

  virtual void tearDown()
  {
    this->clear();
    delete _ctx;
    delete this->_example;
  }

};

class ArrheniusRateVexCLDoubleTest : public ArrheniusRateVectorTestBase<vex::vector<double> >
{
public:

  CPPUNIT_TEST_SUITE( ArrheniusRateVexCLDoubleTest );

  CPPUNIT_TEST( test_standard_rate );
  CPPUNIT_TEST( test_standard_deriv );
  CPPUNIT_TEST( test_standard_rate_and_deriv );

  CPPUNIT_TEST_SUITE_END();

  private:

  vex::Context* _ctx;

public:

  virtual void setUp()
  {
    this->init();
    _ctx = new vex::Context(vex::Filter::DoublePrecision);
    this->_example = new vex::vector<double>(*_ctx, 2*ANTIOCH_N_TUPLES);
  }

  virtual void tearDown()
  {
    this->clear();
    delete _ctx;
    delete this->_example;
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( ArrheniusRateVexCLFloatTest );
CPPUNIT_TEST_SUITE_REGISTRATION( ArrheniusRateVexCLDoubleTest );

#endif // ANTIOCH_HAVE_VEXCL

#endif // ANTIOCH_HAVE_CPPUNIT
