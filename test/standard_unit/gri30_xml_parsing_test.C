#include "antioch_config.h"

#ifdef ANTIOCH_HAVE_CPPUNIT

// CppUnit
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

// C++
#include <limits>

// Antioch
#include "nasa7_thermo_test_base.h"
#include "antioch/chemical_mixture.h"
#include "antioch/nasa7_curve_fit.h"
#include "antioch/nasa_mixture.h"
#include "antioch/nasa_mixture_parsing.h"
#include "antioch/xml_parser.h"
#include "antioch/reaction_set.h"
#include "antioch/read_reaction_set_data.h"

namespace AntiochTesting
{
  template<typename Scalar>
  class GRI30XMLParsingTest : public NASA7ThermoTestBase<Scalar>,
                              public CppUnit::TestCase
  {
  public:

    virtual void setUp()
    {
      this->init();
    }

 
    void test_gri30_xml()
    {
      std::string thermo_filename = std::string(ANTIOCH_SHARE_XML_INPUT_FILES_SOURCE_PATH)+"gri30.xml";      
      Antioch::XMLParser<Scalar> xml_parser(thermo_filename,false);
      std::vector<std::string> species_str_list = xml_parser.species_list();
      
      Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );
      Antioch::NASAThermoMixture<Scalar, Antioch::NASA7CurveFit<Scalar> > nasa_mixture( chem_mixture );

      //xml_parser.read_thermodynamic_data(nasa_mixture);
      Antioch::read_nasa_mixture_data( nasa_mixture, thermo_filename, Antioch::XML );

      Antioch::ReactionSet<Scalar> reaction_set( chem_mixture );
      Antioch::read_reaction_set_data_xml<Scalar>(thermo_filename, true, reaction_set);

      this->check_curve_fits(nasa_mixture);
     }
    
    //could add to check all element coefficients 
    void check_curve_fits( const Antioch::NASAThermoMixture<Scalar, Antioch::NASA7CurveFit<Scalar> >& nasa_mixture)
    {
      const Antioch::NASA7CurveFit<Scalar>& H2_curve_fit =  nasa_mixture.curve_fit(0);
      const Antioch::NASA7CurveFit<Scalar>& N2_curve_fit =  nasa_mixture.curve_fit(47);

      CPPUNIT_ASSERT_EQUAL( (unsigned int)2, H2_curve_fit.n_intervals() );
      CPPUNIT_ASSERT_EQUAL( (unsigned int)2, N2_curve_fit.n_intervals() );

      // Check H2 coefficients
      this->check_coefficients( H2_curve_fit.coefficients(0), this->_H2_coeffs_200_1000 );
      this->check_coefficients( H2_curve_fit.coefficients(1), this->_H2_coeffs_1000_3500 );

      // Check N2 coefficients
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


  

#define DEFINE_GRI30XMLPARSING_SCALAR_TEST(Classname,BaseClass,Scalar)  \
  class Classname : public BaseClass<Scalar>                            \
  {                                                                     \
  public:                                                               \
    CPPUNIT_TEST_SUITE( Classname );                                    \
    CPPUNIT_TEST(test_gri30_xml);                                       \
    CPPUNIT_TEST_SUITE_END();                                           \
  }

  DEFINE_GRI30XMLPARSING_SCALAR_TEST(GRI30XMLParsingTestFloat,GRI30XMLParsingTest,float);
  DEFINE_GRI30XMLPARSING_SCALAR_TEST(GRI30XMLParsingTestDouble,GRI30XMLParsingTest,double);
  DEFINE_GRI30XMLPARSING_SCALAR_TEST(GRI30XMLParsingTestLongDouble,GRI30XMLParsingTest,long double);


  CPPUNIT_TEST_SUITE_REGISTRATION( GRI30XMLParsingTestFloat );
  CPPUNIT_TEST_SUITE_REGISTRATION( GRI30XMLParsingTestDouble );
  CPPUNIT_TEST_SUITE_REGISTRATION( GRI30XMLParsingTestLongDouble );

} // end namespace AntiochTesting



#endif // ANTIOCH_HAVE_CPPUNIT
