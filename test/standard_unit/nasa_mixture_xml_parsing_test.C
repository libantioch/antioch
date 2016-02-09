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
#include "nasa7_thermo_test_base.h"
#include "nasa9_thermo_test_base.h"
#include "antioch/chemical_mixture.h"
#include "antioch/nasa7_curve_fit.h"
#include "antioch/nasa9_curve_fit.h"
#include "antioch/nasa_mixture.h"
#include "antioch/nasa_mixture_parsing.h"

namespace AntiochTesting
{
  template<typename Scalar>
  class NASA9XMLParsingTest : public NASA9ThermoTestBase<Scalar>
  {
  public:

    void test_supplied_species()
    {
      std::vector<std::string> species_str_list(2);
      species_str_list[0] = "N2";
      species_str_list[1] = "NO2";

      Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );
      Antioch::NASAThermoMixture<Scalar, Antioch::NASA9CurveFit<Scalar> > nasa_mixture( chem_mixture );

      std::string filename = std::string(ANTIOCH_SHARE_XML_INPUT_FILES_SOURCE_PATH)+"nasa9_thermo.xml";

      Antioch::read_nasa_mixture_data( nasa_mixture, filename, Antioch::XML );

      const Antioch::NASA9CurveFit<Scalar>& N2_curve_fit =  nasa_mixture.curve_fit(0);
      const Antioch::NASA9CurveFit<Scalar>& NO2_curve_fit =  nasa_mixture.curve_fit(1);

      CPPUNIT_ASSERT_EQUAL( (unsigned int)5, N2_curve_fit.n_intervals() );
      CPPUNIT_ASSERT_EQUAL( (unsigned int)4, NO2_curve_fit.n_intervals() );

      // Check N2 coefficients
      this->check_coefficients( N2_curve_fit.coefficients(0), this->_N2_coeffs_0_200 );
      this->check_coefficients( N2_curve_fit.coefficients(1), this->_N2_coeffs_200_1000 );
      this->check_coefficients( N2_curve_fit.coefficients(2), this->_N2_coeffs_1000_6000 );
      this->check_coefficients( N2_curve_fit.coefficients(3), this->_N2_coeffs_6000_20000 );
      this->check_coefficients( N2_curve_fit.coefficients(4), this->_N2_coeffs_20000_99999 );

      // Check NO2 coefficients
      this->check_coefficients( NO2_curve_fit.coefficients(0), this->_NO2_coeffs_0_200 );
      this->check_coefficients( NO2_curve_fit.coefficients(1), this->_NO2_coeffs_200_1000 );
      this->check_coefficients( NO2_curve_fit.coefficients(2), this->_NO2_coeffs_1000_6000 );
      this->check_coefficients( NO2_curve_fit.coefficients(3), this->_NO2_coeffs_6000_99999 );
    }

    void check_coefficients( const Scalar* parsed_coeffs, std::vector<Scalar>& exact_coeffs )
    {
      for( unsigned int c = 0; c < exact_coeffs.size(); c++ )
        CPPUNIT_ASSERT_DOUBLES_EQUAL( exact_coeffs[c], parsed_coeffs[c], this->tol() );
    }

    Scalar tol()
    { return std::numeric_limits<Scalar>::epsilon() * 10; }

  };

  template<typename Scalar>
  class NASA7XMLParsingTest : public NASA7ThermoTestBase<Scalar>
  {
  public:

    void test_supplied_species()
    {
      std::vector<std::string> species_str_list(2);
      species_str_list[0] = "H2";
      species_str_list[1] = "N2";

      Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );
      Antioch::NASAThermoMixture<Scalar, Antioch::NASA7CurveFit<Scalar> > nasa_mixture( chem_mixture );

      std::string filename = std::string(ANTIOCH_TESTING_INPUT_FILES_PATH)+"nasa7_thermo_test.xml";

      Antioch::read_nasa_mixture_data( nasa_mixture, filename, Antioch::XML );

      const Antioch::NASA7CurveFit<Scalar>& H2_curve_fit =  nasa_mixture.curve_fit(0);
      const Antioch::NASA7CurveFit<Scalar>& N2_curve_fit =  nasa_mixture.curve_fit(1);

      CPPUNIT_ASSERT_EQUAL( (unsigned int)2, H2_curve_fit.n_intervals() );
      CPPUNIT_ASSERT_EQUAL( (unsigned int)2, N2_curve_fit.n_intervals() );

      // Check N2 coefficients
      this->check_coefficients( H2_curve_fit.coefficients(0), this->_H2_coeffs_200_1000 );
      this->check_coefficients( H2_curve_fit.coefficients(1), this->_H2_coeffs_1000_3500 );

      // Check NO2 coefficients
      this->check_coefficients( N2_curve_fit.coefficients(0), this->_N2_coeffs_300_1000 );
      this->check_coefficients( N2_curve_fit.coefficients(1), this->_N2_coeffs_1000_5000 );
     }

    void check_coefficients( const Scalar* parsed_coeffs, std::vector<Scalar>& exact_coeffs )
    {
      for( unsigned int c = 0; c < exact_coeffs.size(); c++ )
        CPPUNIT_ASSERT_DOUBLES_EQUAL( exact_coeffs[c], parsed_coeffs[c], this->tol() );
    }

    Scalar tol()
    { return std::numeric_limits<Scalar>::epsilon() * 10; }

  };

  class NASA9XMLParsingTestFloat : public NASA9XMLParsingTest<float>
  {
  public:
    CPPUNIT_TEST_SUITE( NASA9XMLParsingTestFloat );

    CPPUNIT_TEST(test_supplied_species);

    CPPUNIT_TEST_SUITE_END();

  };

  class NASA9XMLParsingTestDouble : public NASA9XMLParsingTest<double>
  {
  public:
    CPPUNIT_TEST_SUITE( NASA9XMLParsingTestDouble );

    CPPUNIT_TEST(test_supplied_species);

    CPPUNIT_TEST_SUITE_END();

  };

  class NASA9XMLParsingTestLongDouble : public NASA9XMLParsingTest<long double>
  {
  public:
    CPPUNIT_TEST_SUITE( NASA9XMLParsingTestLongDouble );

    CPPUNIT_TEST(test_supplied_species);

    CPPUNIT_TEST_SUITE_END();

  };

  class NASA7XMLParsingTestFloat : public NASA7XMLParsingTest<float>
  {
  public:
    CPPUNIT_TEST_SUITE( NASA7XMLParsingTestFloat );

    CPPUNIT_TEST(test_supplied_species);

    CPPUNIT_TEST_SUITE_END();

  };

  class NASA7XMLParsingTestDouble : public NASA7XMLParsingTest<double>
  {
  public:
    CPPUNIT_TEST_SUITE( NASA7XMLParsingTestDouble );

    CPPUNIT_TEST(test_supplied_species);

    CPPUNIT_TEST_SUITE_END();

  };

  class NASA7XMLParsingTestLongDouble : public NASA7XMLParsingTest<long double>
  {
  public:
    CPPUNIT_TEST_SUITE( NASA7XMLParsingTestLongDouble );

    CPPUNIT_TEST(test_supplied_species);

    CPPUNIT_TEST_SUITE_END();

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( NASA9XMLParsingTestFloat );
  CPPUNIT_TEST_SUITE_REGISTRATION( NASA9XMLParsingTestDouble );
  CPPUNIT_TEST_SUITE_REGISTRATION( NASA9XMLParsingTestLongDouble );
  CPPUNIT_TEST_SUITE_REGISTRATION( NASA7XMLParsingTestFloat );
  CPPUNIT_TEST_SUITE_REGISTRATION( NASA7XMLParsingTestDouble );
  CPPUNIT_TEST_SUITE_REGISTRATION( NASA7XMLParsingTestLongDouble );

} // end namespace AntiochTesting

#endif // ANTIOCH_HAVE_CPPUNIT
