#include "antioch_config.h"

#ifdef ANTIOCH_HAVE_CPPUNIT

// CppUnit
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

// C++
#include <limits>

// Antioch
#include "antioch/xml_parser.h"
#include "antioch/reaction_set.h"
#include "antioch/read_reaction_set_data.h"

namespace AntiochTesting
{
  template<typename Scalar>
  class XMLMultiSectionParsingTestBase : public CppUnit::TestCase
  {
  public:

    void test_ozone()
    {
      std::string filename = std::string(ANTIOCH_TESTING_INPUT_FILES_PATH)+"multiple_secs.xml";
      const std::string mixture("ozone");

      Antioch::XMLParser<Scalar> parser(filename,mixture,false);

      std::vector<std::string> species_names = parser.species_list();
      this->check_ozone_species(species_names);

      Antioch::ChemicalMixture<Scalar> chem_mixture( species_names );
      Antioch::ReactionSet<Scalar> reaction_set( chem_mixture );

      Antioch::read_reaction_set_data<Scalar>(false,reaction_set,&parser);
      this->check_ozone_reactions(reaction_set);
    }

    void test_air5sp()
    {
      std::string filename = std::string(ANTIOCH_TESTING_INPUT_FILES_PATH)+"multiple_secs.xml";
      const std::string mixture("air5sp");

      Antioch::XMLParser<Scalar> parser(filename,mixture,false);

      std::vector<std::string> species_names = parser.species_list();
      this->check_air5sp_species(species_names);

      Antioch::ChemicalMixture<Scalar> chem_mixture( species_names );
      Antioch::ReactionSet<Scalar> reaction_set( chem_mixture );

      Antioch::read_reaction_set_data<Scalar>(false,reaction_set,&parser);
      this->check_air5sp_reactions(reaction_set);
    }

  protected:

    void check_ozone_species( const std::vector<std::string> & species_names )
    {
      CPPUNIT_ASSERT_EQUAL(3,(int)species_names.size());
      CPPUNIT_ASSERT_EQUAL(std::string("O"),species_names[0]);
      CPPUNIT_ASSERT_EQUAL(std::string("O2"),species_names[1]);
      CPPUNIT_ASSERT_EQUAL(std::string("O3"),species_names[2]);
    }

    void check_air5sp_species( const std::vector<std::string> & species_names )
    {
      CPPUNIT_ASSERT_EQUAL(5,(int)species_names.size());
      CPPUNIT_ASSERT_EQUAL(std::string("N2"),species_names[0]);
      CPPUNIT_ASSERT_EQUAL(std::string("O2"),species_names[1]);
      CPPUNIT_ASSERT_EQUAL(std::string("NO"),species_names[2]);
      CPPUNIT_ASSERT_EQUAL(std::string("N"),species_names[3]);
      CPPUNIT_ASSERT_EQUAL(std::string("O"),species_names[4]);
    }

    void check_ozone_reactions( const Antioch::ReactionSet<Scalar> & reaction_set )
    {
      CPPUNIT_ASSERT_EQUAL(3,(int)reaction_set.n_reactions());

      const Antioch::ChemicalMixture<Scalar> &  mixture = reaction_set.chemical_mixture();

      // Check first reaction
      {
        const Antioch::Reaction<Scalar> & modarr_threebody = reaction_set.reaction(0);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::THREE_BODY, modarr_threebody.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::KOOIJ, modarr_threebody.kinetics_model() );
        CPPUNIT_ASSERT(modarr_threebody.reversible());

        // Check threebody efficiencies
        for( unsigned int s = 0; s < reaction_set.n_species(); s++ )
          {
            std::string species_name = mixture.chemical_species()[s]->species();
            Scalar efficiency = modarr_threebody.get_efficiency(s);

            if( species_name == std::string("O3") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0.92, efficiency, this->tol());
            else if( species_name == std::string("O2") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0.94, efficiency, this->tol());
            else if( species_name == std::string("O") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(1.13, efficiency, this->tol());
            else
              // Default should be 1.0
              CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, efficiency, this->tol());
          }
      }


      // Check second reaction
      {
        const Antioch::Reaction<Scalar> & arr_threebody = reaction_set.reaction(1);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::THREE_BODY, arr_threebody.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::ARRHENIUS, arr_threebody.kinetics_model() );
        CPPUNIT_ASSERT(arr_threebody.reversible());

        // Check threebody efficiencies
        for( unsigned int s = 0; s < reaction_set.n_species(); s++ )
          {
            std::string species_name = mixture.chemical_species()[s]->species();
            Scalar efficiency = arr_threebody.get_efficiency(s);

            if( species_name == std::string("O3") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0.92, efficiency, this->tol());
            else if( species_name == std::string("O2") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0.94, efficiency, this->tol());
            else if( species_name == std::string("O") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(1.13, efficiency, this->tol());
            else
              // Default should be 1.0
              CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, efficiency, this->tol());
          }
      }

      // Check third reaction
      {
        const Antioch::Reaction<Scalar> & arr_rev = reaction_set.reaction(2);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::ELEMENTARY ,arr_rev.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::ARRHENIUS, arr_rev.kinetics_model() );
        CPPUNIT_ASSERT(arr_rev.reversible());
      }
    }


    void check_air5sp_reactions( const Antioch::ReactionSet<Scalar> & reaction_set )
    {
      const Antioch::ChemicalMixture<Scalar> &  mixture = reaction_set.chemical_mixture();

      // This file actually has a lot more reactions, so we're additionally
      // testing that we're stripping off reactions that include species
      // that are not in the mixture that formed the reaction set.
      CPPUNIT_ASSERT_EQUAL(5,(int)reaction_set.n_reactions());

      // Check first reaction
      {
        const Antioch::Reaction<Scalar> & modarr_threebody = reaction_set.reaction(0);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::THREE_BODY, modarr_threebody.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::KOOIJ, modarr_threebody.kinetics_model() );
        CPPUNIT_ASSERT(modarr_threebody.reversible());

        // Check threebody efficiencies
        for( unsigned int s = 0; s < reaction_set.n_species(); s++ )
          {
            std::string species_name = mixture.chemical_species()[s]->species();
            Scalar efficiency = modarr_threebody.get_efficiency(s);

            if( species_name == std::string("N") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(4.2857, efficiency, this->tol());
            else if( species_name == std::string("O") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(4.2857, efficiency, this->tol());
            else
              // Default should be 1.0
              CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, efficiency, this->tol());
          }
      }

      // Check second reaction
      {
        const Antioch::Reaction<Scalar> & modarr_threebody = reaction_set.reaction(1);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::THREE_BODY, modarr_threebody.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::KOOIJ, modarr_threebody.kinetics_model() );
        CPPUNIT_ASSERT(modarr_threebody.reversible());

        // Check threebody efficiencies
        for( unsigned int s = 0; s < reaction_set.n_species(); s++ )
          {
            std::string species_name = mixture.chemical_species()[s]->species();
            Scalar efficiency = modarr_threebody.get_efficiency(s);

            if( species_name == std::string("N") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, efficiency, this->tol());
            else if( species_name == std::string("O") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, efficiency, this->tol());
            else
              // Default should be 1.0
              CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, efficiency, this->tol());
          }
      }

      // Check third reaction
      {
        const Antioch::Reaction<Scalar> & arr_threebody = reaction_set.reaction(2);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::THREE_BODY, arr_threebody.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::ARRHENIUS, arr_threebody.kinetics_model() );
        CPPUNIT_ASSERT(arr_threebody.reversible());

        // Check threebody efficiencies
        for( unsigned int s = 0; s < reaction_set.n_species(); s++ )
          {
            std::string species_name = mixture.chemical_species()[s]->species();
            Scalar efficiency = arr_threebody.get_efficiency(s);

            if( species_name == std::string("N") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(22.0, efficiency, this->tol());
            else if( species_name == std::string("O") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(22.0, efficiency, this->tol());
            else if( species_name == std::string("NO") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(22.0, efficiency, this->tol());
            else
              // Default should be 1.0
              CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, efficiency, this->tol());
          }
      }


      // Check fourth reaction
      {
        const Antioch::Reaction<Scalar> & modarr_rev = reaction_set.reaction(3);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::ELEMENTARY ,modarr_rev.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::KOOIJ, modarr_rev.kinetics_model() );
        CPPUNIT_ASSERT(modarr_rev.reversible());
      }

      // Check fifth reaction
      {
        const Antioch::Reaction<Scalar> & arr_rev = reaction_set.reaction(4);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::ELEMENTARY ,arr_rev.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::ARRHENIUS, arr_rev.kinetics_model() );
        CPPUNIT_ASSERT(arr_rev.reversible());
      }

    }

    Scalar tol()
    { return std::numeric_limits<Scalar>::epsilon() * 10; }

  };

  #define DEFINE_XMLPARSING_SCALAR_TEST(Classname,BaseClass,Scalar)  \
  class Classname : public BaseClass<Scalar>                            \
  {                                                                     \
  public:                                                               \
    CPPUNIT_TEST_SUITE( Classname );                                    \
    CPPUNIT_TEST(test_ozone);                                           \
    CPPUNIT_TEST(test_air5sp);                                          \
    CPPUNIT_TEST_SUITE_END();                                           \
  }

  DEFINE_XMLPARSING_SCALAR_TEST(XMLMultiSectionParsingTestFloat,XMLMultiSectionParsingTestBase,float);
  DEFINE_XMLPARSING_SCALAR_TEST(XMLMultiSectionParsingTestDouble,XMLMultiSectionParsingTestBase,double);
  DEFINE_XMLPARSING_SCALAR_TEST(XMLMultiSectionParsingTestLongDouble,XMLMultiSectionParsingTestBase,long double);


  CPPUNIT_TEST_SUITE_REGISTRATION( XMLMultiSectionParsingTestFloat );
  CPPUNIT_TEST_SUITE_REGISTRATION( XMLMultiSectionParsingTestDouble );
  CPPUNIT_TEST_SUITE_REGISTRATION( XMLMultiSectionParsingTestLongDouble );

} // end namespace AntiochTesting



#endif // ANTIOCH_HAVE_CPPUNIT
