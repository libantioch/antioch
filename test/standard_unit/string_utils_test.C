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

// Antioch
#include "antioch/string_utils.h"

namespace AntiochTesting
{
  class StringUtilitiesTest : public CppUnit::TestCase
  {
  public:

    void setUp(){}

    void tearDown(){}

    CPPUNIT_TEST_SUITE( StringUtilitiesTest );

    CPPUNIT_TEST(test_split_string);

    CPPUNIT_TEST_SUITE_END();

  private:

    void test_split_string()
    {
      {
        std::string str("N->N2");
        std::vector<std::string> str_split;
        Antioch::split_string( str, "->", str_split);

        CPPUNIT_ASSERT_EQUAL(2,(int)str_split.size());
        CPPUNIT_ASSERT_EQUAL(str_split[0],std::string("N"));
        CPPUNIT_ASSERT_EQUAL(str_split[1],std::string("N2"));
      }

      {
        std::string str("N+C(s)->CN");
        std::vector<std::string> str_split;
        Antioch::split_string( str, "->", str_split);

        CPPUNIT_ASSERT_EQUAL(2,(int)str_split.size());
        CPPUNIT_ASSERT_EQUAL(str_split[0],std::string("N+C(s)"));
        CPPUNIT_ASSERT_EQUAL(str_split[1],std::string("CN"));
      }

      {
        std::string str("u:v:w:T:p:w_N:w_N2:p0");
        std::vector<std::string> str_split;
        Antioch::split_string( str, ":", str_split);

        CPPUNIT_ASSERT_EQUAL(8,(int)str_split.size());
        CPPUNIT_ASSERT_EQUAL(str_split[0],std::string("u"));
        CPPUNIT_ASSERT_EQUAL(str_split[1],std::string("v"));
        CPPUNIT_ASSERT_EQUAL(str_split[2],std::string("w"));
        CPPUNIT_ASSERT_EQUAL(str_split[3],std::string("T"));
        CPPUNIT_ASSERT_EQUAL(str_split[4],std::string("p"));
        CPPUNIT_ASSERT_EQUAL(str_split[5],std::string("w_N"));
        CPPUNIT_ASSERT_EQUAL(str_split[6],std::string("w_N2"));
        CPPUNIT_ASSERT_EQUAL(str_split[7],std::string("p0"));
      }

      {
        std::string str("u v w T p w_N w_N2 p0");
        std::vector<std::string> str_split;
        Antioch::split_string( str, " ", str_split);

        CPPUNIT_ASSERT_EQUAL(8,(int)str_split.size());
        CPPUNIT_ASSERT_EQUAL(str_split[0],std::string("u"));
        CPPUNIT_ASSERT_EQUAL(str_split[1],std::string("v"));
        CPPUNIT_ASSERT_EQUAL(str_split[2],std::string("w"));
        CPPUNIT_ASSERT_EQUAL(str_split[3],std::string("T"));
        CPPUNIT_ASSERT_EQUAL(str_split[4],std::string("p"));
        CPPUNIT_ASSERT_EQUAL(str_split[5],std::string("w_N"));
        CPPUNIT_ASSERT_EQUAL(str_split[6],std::string("w_N2"));
        CPPUNIT_ASSERT_EQUAL(str_split[7],std::string("p0"));
      }

      {
        std::string str("2.25853318e+03,    -1.57454401e+00,    2.50363759e+00,    -5.20253954e-06, ");
        str += "   4.51647956e-09,  -2.18115104e-12,   4.49430845e-16, 2.16895191e+05, 4.34577527e+00";

        std::vector<std::string> str_split;
        Antioch::split_string( str, " ,", str_split);

        CPPUNIT_ASSERT_EQUAL(9,(int)str_split.size());
        CPPUNIT_ASSERT_EQUAL(str_split[0],std::string("2.25853318e+03"));
        CPPUNIT_ASSERT_EQUAL(str_split[1],std::string("-1.57454401e+00"));
        CPPUNIT_ASSERT_EQUAL(str_split[2],std::string("2.50363759e+00"));
        CPPUNIT_ASSERT_EQUAL(str_split[3],std::string("-5.20253954e-06"));
        CPPUNIT_ASSERT_EQUAL(str_split[4],std::string("4.51647956e-09"));
        CPPUNIT_ASSERT_EQUAL(str_split[5],std::string("-2.18115104e-12"));
        CPPUNIT_ASSERT_EQUAL(str_split[6],std::string("4.49430845e-16"));
        CPPUNIT_ASSERT_EQUAL(str_split[7],std::string("2.16895191e+05"));
        CPPUNIT_ASSERT_EQUAL(str_split[8],std::string("4.34577527e+00"));
      }
    }
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( StringUtilitiesTest );

} // end namespace AntiochTesting

#endif // ANTIOCH_HAVE_CPPUNIT
