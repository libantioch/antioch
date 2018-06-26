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
      this->init_exact_species();
    }

 
    void test_gri30_xml()
    {
      std::string thermo_filename = std::string(ANTIOCH_SHARE_XML_INPUT_FILES_SOURCE_PATH)+"gri30.xml";
      const std::string phase("gri30_mix");

      Antioch::XMLParser<Scalar> xml_parser(thermo_filename,phase,false);
      std::vector<std::string> species_str_list = xml_parser.species_list();

      this->check_species_list(species_str_list);

      Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );
      Antioch::NASAThermoMixture<Scalar, Antioch::NASA7CurveFit<Scalar> > nasa_mixture( chem_mixture );

      //xml_parser.read_thermodynamic_data(nasa_mixture);
      Antioch::read_nasa_mixture_data( nasa_mixture, thermo_filename, Antioch::XML );
      this->check_curve_fits(nasa_mixture);

      Antioch::ReactionSet<Scalar> reaction_set( chem_mixture );
      Antioch::read_reaction_set_data_xml<Scalar>(thermo_filename, true, reaction_set);
      this->check_reaction_set(reaction_set);
     }

    void check_species_list(const std::vector<std::string> & species_list)
    {
      CPPUNIT_ASSERT_EQUAL(_species_exact.size(), species_list.size());

      // Check that these particular species are where they are supposed to be
      // since we check them in the thermo part as well
      CPPUNIT_ASSERT_EQUAL( std::string("H2"), species_list[_H2_species_id] );
      CPPUNIT_ASSERT_EQUAL( std::string("N2"), species_list[_N2_species_id] );
      CPPUNIT_ASSERT_EQUAL( std::string("HCNO"), species_list[_HCNO_species_id] );

      for( unsigned int s = 0; s < _species_exact.size(); s++ )
        CPPUNIT_ASSERT_EQUAL( _species_exact[s], species_list[s] );
    }

    //could add to check all element coefficients
    void check_curve_fits( const Antioch::NASAThermoMixture<Scalar, Antioch::NASA7CurveFit<Scalar> >& nasa_mixture)
    {
      const Antioch::NASA7CurveFit<Scalar> & H2_curve_fit =
        nasa_mixture.curve_fit(_H2_species_id);

      const Antioch::NASA7CurveFit<Scalar> & N2_curve_fit =
        nasa_mixture.curve_fit(_N2_species_id);

      const Antioch::NASA7CurveFit<Scalar> & HCNO_curve_fit =
        nasa_mixture.curve_fit(_HCNO_species_id);

      CPPUNIT_ASSERT_EQUAL( (unsigned int)2, H2_curve_fit.n_intervals() );
      CPPUNIT_ASSERT_EQUAL( (unsigned int)2, N2_curve_fit.n_intervals() );
      CPPUNIT_ASSERT_EQUAL( (unsigned int)2, HCNO_curve_fit.n_intervals() );

      // Check H2 coefficients
      this->check_coefficients( H2_curve_fit.coefficients(0), this->_H2_coeffs_200_1000 );
      this->check_coefficients( H2_curve_fit.coefficients(1), this->_H2_coeffs_1000_3500 );

      // Check N2 coefficients
      this->check_coefficients( N2_curve_fit.coefficients(0), this->_N2_coeffs_300_1000 );
      this->check_coefficients( N2_curve_fit.coefficients(1), this->_N2_coeffs_1000_5000 );

      // Check the bounds on HCNO curve fit since they're non-standard
      // We don't give access to the temperature bounds, so we just check
      // that the interval changes when T crosses the critical value
      CPPUNIT_ASSERT_EQUAL( 0, (int)HCNO_curve_fit.interval(1381.0) );
      CPPUNIT_ASSERT_EQUAL( 1, (int)HCNO_curve_fit.interval(1383.0) );
    }

    void check_coefficients( const Scalar* parsed_coeffs, std::vector<Scalar>& exact_coeffs )
    {
      for( unsigned int c = 0; c < exact_coeffs.size(); c++ )
        CPPUNIT_ASSERT_DOUBLES_EQUAL( exact_coeffs[c], parsed_coeffs[c], this->tol() );
    }

    Scalar tol()
    { return std::numeric_limits<Scalar>::epsilon() * 10; }

    void check_reaction_set( const Antioch::ReactionSet<Scalar> & reaction_set )
    {
      CPPUNIT_ASSERT_EQUAL( (int)_species_exact.size(), (int)reaction_set.n_species() );
      CPPUNIT_ASSERT_EQUAL( 325, (int)reaction_set.n_reactions() );
    }

  protected:

    std::vector<std::string> _species_exact;
    unsigned int _H2_species_id;
    unsigned int _N2_species_id;
    unsigned int _HCNO_species_id;

    void init_exact_species()
    {
      _H2_species_id = 0;
      _N2_species_id = 47;
      _HCNO_species_id = 43;

      _species_exact.reserve(53);
      _species_exact.push_back("H2");
      _species_exact.push_back("H");
      _species_exact.push_back("O");
      _species_exact.push_back("O2");
      _species_exact.push_back("OH");
      _species_exact.push_back("H2O");
      _species_exact.push_back("HO2");
      _species_exact.push_back("H2O2");
      _species_exact.push_back("C");
      _species_exact.push_back("CH");
      _species_exact.push_back("CH2");
      _species_exact.push_back("CH2(S)");
      _species_exact.push_back("CH3");
      _species_exact.push_back("CH4");
      _species_exact.push_back("CO");
      _species_exact.push_back("CO2");
      _species_exact.push_back("HCO");
      _species_exact.push_back("CH2O");
      _species_exact.push_back("CH2OH");
      _species_exact.push_back("CH3O");
      _species_exact.push_back("CH3OH");
      _species_exact.push_back("C2H");
      _species_exact.push_back("C2H2");
      _species_exact.push_back("C2H3");
      _species_exact.push_back("C2H4");
      _species_exact.push_back("C2H5");
      _species_exact.push_back("C2H6");
      _species_exact.push_back("HCCO");
      _species_exact.push_back("CH2CO");
      _species_exact.push_back("HCCOH");
      _species_exact.push_back("N");
      _species_exact.push_back("NH");
      _species_exact.push_back("NH2");
      _species_exact.push_back("NH3");
      _species_exact.push_back("NNH");
      _species_exact.push_back("NO");
      _species_exact.push_back("NO2");
      _species_exact.push_back("N2O");
      _species_exact.push_back("HNO");
      _species_exact.push_back("CN");
      _species_exact.push_back("HCN");
      _species_exact.push_back("H2CN");
      _species_exact.push_back("HCNN");
      _species_exact.push_back("HCNO");
      _species_exact.push_back("HOCN");
      _species_exact.push_back("HNCO");
      _species_exact.push_back("NCO");
      _species_exact.push_back("N2");
      _species_exact.push_back("AR");
      _species_exact.push_back("C3H7");
      _species_exact.push_back("C3H8");
      _species_exact.push_back("CH2CHO");
      _species_exact.push_back("CH3CHO");
    }

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
