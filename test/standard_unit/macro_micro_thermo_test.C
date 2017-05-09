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

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

// C++
#include <limits>

// Anitoch testing
#include "testing_utils.h"
#include "nasa_test_helper.h"

// Antioch
#include "antioch/vector_utils_decl.h"

#include "antioch/chemical_mixture.h"
#include "antioch/stat_mech_thermo.h"
#include "antioch/ideal_gas_micro_thermo.h"
#include "antioch/nasa_mixture.h"
#include "antioch/nasa_evaluator.h"
#include "antioch/nasa7_curve_fit.h"
#include "antioch/nasa9_curve_fit.h"
#include "antioch/nasa_mixture_parsing.h"

#include "antioch/vector_utils.h"

namespace AntiochTesting
{
  template<typename Scalar>
  class MacroMicroThermoTestBase
  {
  protected:

    void init()
    {
      _n_species = 5;

      _species_name_list.reserve(_n_species);
      _species_name_list.push_back( "N2" );
      _species_name_list.push_back( "O2" );
      _species_name_list.push_back( "N" );
      _species_name_list.push_back( "O" );
      _species_name_list.push_back( "NO" );

      _molar_mass.resize(5);
      _molar_mass[2]  = 14.008e-3L; //in SI kg/mol
      _molar_mass[3]  = 16e-3L;     //in SI kg/mol
      _molar_mass[0] = 2.L * _molar_mass[2]; //in SI kg/mol
      _molar_mass[1] = 2.L * _molar_mass[3]; //in SI kg/mol
      _molar_mass[4] = _molar_mass[3] + _molar_mass[2]; //in SI kg/mol

      _gas_consts.resize(5);
      _gas_consts[0] = Antioch::Constants::R_universal<Scalar>() / _molar_mass[0];
      _gas_consts[1] = Antioch::Constants::R_universal<Scalar>() / _molar_mass[1];
      _gas_consts[2] = Antioch::Constants::R_universal<Scalar>() / _molar_mass[2];
      _gas_consts[3] = Antioch::Constants::R_universal<Scalar>() / _molar_mass[3];
      _gas_consts[4] = Antioch::Constants::R_universal<Scalar>() / _molar_mass[4];

      _n_tr_dofs.resize(5);
      _n_tr_dofs[0] = 2.5L;
      _n_tr_dofs[1] = 2.5L;
      _n_tr_dofs[2] = 1.5L;
      _n_tr_dofs[3] = 1.5L;
      _n_tr_dofs[4] = 2.5L;

      _chem_mixture = new Antioch::ChemicalMixture<Scalar>( _species_name_list );

      // Mass fractions
      _mass_fractions.resize( _n_species, 0.2 );
      _mass_fractions[0] = 0.5L;
      _mass_fractions[1] = 0.2L;
      _mass_fractions[2] = 0.1L;
      _mass_fractions[3] = 0.1L;
      _mass_fractions[4] = 0.1L;
    }

    void clear()
    {
      delete _chem_mixture;
    }

    unsigned int _n_species;

    std::vector<std::string> _species_name_list;

    std::vector<Scalar> _molar_mass;
    std::vector<Scalar> _gas_consts;
    std::vector<Scalar> _n_tr_dofs;

    std::vector<Scalar> _mass_fractions;

    Antioch::ChemicalMixture<Scalar> * _chem_mixture;
  };

  template<typename Scalar>
  class StatMechThermoTestBase : public CppUnit::TestCase,
                                 public TestingUtilities<Scalar>,
                                 public MacroMicroThermoTestBase<Scalar>
  {
  public:

    void test_cv_trans()
    {
      for( unsigned int s = 0; s < this->_n_species; s++ )
        this->test_scalar_rel( Scalar(1.5)*this->_gas_consts[s],
                               _thermo->cv_trans(s),
                               std::numeric_limits<Scalar>::epsilon()*10 );
    }

    void test_cv_trans_over_R()
    {
      for( unsigned int s = 0; s < this->_n_species; s++ )
        this->test_scalar_rel( Scalar(1.5),
                               _thermo->cv_trans_over_R(s),
                               std::numeric_limits<Scalar>::epsilon()*10 );
    }

    void test_cv_tr()
    {
      for( unsigned int s = 0; s < this->_n_species; s++ )
        this->test_scalar_rel( this->_n_tr_dofs[s]*this->_gas_consts[s],
                               _thermo->cv_tr(s),
                               std::numeric_limits<Scalar>::epsilon()*10 );
    }

    void test_cv_tr_over_R()
    {
      for( unsigned int s = 0; s < this->_n_species; s++ )
        this->test_scalar_rel( this->_n_tr_dofs[s],
                               _thermo->cv_tr_over_R(s),
                               std::numeric_limits<Scalar>::epsilon()*10 );
    }

    void test_cv_rot()
    {
      for( unsigned int s = 0; s < this->_n_species; s++ )
        this->test_scalar_rel( (this->_n_tr_dofs[s]-Scalar(1.5))*this->_gas_consts[s],
                               _thermo->cv_rot(s),
                               std::numeric_limits<Scalar>::epsilon()*10 );
    }

    void test_cv_rot_over_R()
    {
      for( unsigned int s = 0; s < this->_n_species; s++ )
        this->test_scalar_rel( this->_n_tr_dofs[s]-Scalar(1.5),
                               _thermo->cv_rot_over_R(s),
                               std::numeric_limits<Scalar>::epsilon()*10 );
    }

    void test_cv_tr_mix()
    {
      Scalar cv_exact(0.0);

      for( unsigned int s = 0; s < this->_n_species; s++ )
        cv_exact += this->_n_tr_dofs[s]*this->_gas_consts[s]*this->_mass_fractions[s];

      this->test_scalar_rel( cv_exact,
                             _thermo->cv_tr(this->_mass_fractions),
                             std::numeric_limits<Scalar>::epsilon()*10 );
    }

    void setUp()
    {
      this->init();
      _thermo = new Antioch::StatMechThermodynamics<Scalar>( *(this->_chem_mixture) );
    }

    void tearDown()
    {
      this->clear();
      delete _thermo;
    }

  private:

    Antioch::StatMechThermodynamics<Scalar> * _thermo;

  };

#define DEFINE_STATMECHTHERMO_SCALAR_TEST(classname,scalar)     \
  class classname : public StatMechThermoTestBase<scalar>       \
  {                                                             \
  public:                                                       \
    CPPUNIT_TEST_SUITE( classname );                            \
    CPPUNIT_TEST(test_cv_trans);                                \
    CPPUNIT_TEST(test_cv_trans_over_R);                         \
    CPPUNIT_TEST(test_cv_tr);                                   \
    CPPUNIT_TEST(test_cv_tr_over_R);                            \
    CPPUNIT_TEST(test_cv_rot);                                  \
    CPPUNIT_TEST(test_cv_rot_over_R);                           \
    CPPUNIT_TEST(test_cv_tr_mix);                               \
    CPPUNIT_TEST_SUITE_END();                                   \
  }

  DEFINE_STATMECHTHERMO_SCALAR_TEST(StatMechThermoTestFloat,float);
  DEFINE_STATMECHTHERMO_SCALAR_TEST(StatMechThermoTestDouble,double);
  DEFINE_STATMECHTHERMO_SCALAR_TEST(StatMechThermoTestLongDouble,long double);


  CPPUNIT_TEST_SUITE_REGISTRATION( StatMechThermoTestFloat );
  CPPUNIT_TEST_SUITE_REGISTRATION( StatMechThermoTestDouble );
  CPPUNIT_TEST_SUITE_REGISTRATION( StatMechThermoTestLongDouble );




  template<typename Scalar, typename NASACurveFit>
  class IdealGasMicroThermoTestBase : public CppUnit::TestCase,
                                      public TestingUtilities<Scalar>,
                                      public MacroMicroThermoTestBase<Scalar>
  {
  public:

    void test_cv_trans()
    {
      for( unsigned int s = 0; s < this->_n_species; s++ )
        this->test_scalar_rel( Scalar(1.5)*this->_gas_consts[s],
                               _thermo->cv_trans(s),
                               std::numeric_limits<Scalar>::epsilon()*10 );
    }

    void test_cv_trans_over_R()
    {
      for( unsigned int s = 0; s < this->_n_species; s++ )
        this->test_scalar_rel( Scalar(1.5),
                               _thermo->cv_trans_over_R(s),
                               std::numeric_limits<Scalar>::epsilon()*10 );
    }

    void test_cv_tr()
    {
      for( unsigned int s = 0; s < this->_n_species; s++ )
        this->test_scalar_rel( this->_n_tr_dofs[s]*this->_gas_consts[s],
                               _thermo->cv_tr(s),
                               std::numeric_limits<Scalar>::epsilon()*10 );
    }

    void test_cv_tr_over_R()
    {
      for( unsigned int s = 0; s < this->_n_species; s++ )
        this->test_scalar_rel( this->_n_tr_dofs[s],
                               _thermo->cv_tr_over_R(s),
                               std::numeric_limits<Scalar>::epsilon()*10 );
    }

    void test_cv_rot()
    {
      for( unsigned int s = 0; s < this->_n_species; s++ )
        this->test_scalar_rel( (this->_n_tr_dofs[s]-Scalar(1.5))*this->_gas_consts[s],
                               _thermo->cv_rot(s),
                               std::numeric_limits<Scalar>::epsilon()*10 );
    }

    void test_cv_rot_over_R()
    {
      for( unsigned int s = 0; s < this->_n_species; s++ )
        this->test_scalar_rel( this->_n_tr_dofs[s]-Scalar(1.5),
                               _thermo->cv_rot_over_R(s),
                               std::numeric_limits<Scalar>::epsilon()*10 );
    }

    void test_cv_tr_mix()
    {
      Scalar cv_exact(0.0);

      for( unsigned int s = 0; s < this->_n_species; s++ )
        cv_exact += this->_n_tr_dofs[s]*this->_gas_consts[s]*this->_mass_fractions[s];

      this->test_scalar_rel( cv_exact,
                             _thermo->cv_tr(this->_mass_fractions),
                             std::numeric_limits<Scalar>::epsilon()*10 );
    }

    void setUp()
    {
      this->init();
      this->init_nasa_coeffs();

      _nasa_mixture = new Antioch::NASAThermoMixture<Scalar,NASACurveFit>( *(this->_chem_mixture) );

      this->parse_nasa_coeffs( *(this->_nasa_mixture) );


      _nasa_evaluator = new Antioch::NASAEvaluator<Scalar,NASACurveFit>( *(_nasa_mixture) );

      _thermo = new Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<Scalar,NASACurveFit>,Scalar>
        ( *(_nasa_evaluator), *(this->_chem_mixture) );
    }

    void tearDown()
    {
      delete _thermo;
      delete _nasa_evaluator;
      delete _nasa_mixture;
      this->clear();
    }

  protected:

    Antioch::NASAThermoMixture<Scalar,NASACurveFit> * _nasa_mixture;

    Antioch::NASAEvaluator<Scalar,NASACurveFit> * _nasa_evaluator;

    Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<Scalar,NASACurveFit>,Scalar> * _thermo;

    std::map<std::string,std::vector<Scalar> > _nasa_coeffs_lower_interval;
    std::map<std::string,std::vector<Scalar> > _nasa_coeffs_upper_interval;

    virtual void parse_nasa_coeffs( Antioch::NASAThermoMixture<Scalar,NASACurveFit> & nasa_mixture ) =0;

    virtual void init_nasa_coeffs() =0;

  };

  template<typename Scalar>
  class IdealGasMicroThermoTestNASA7Base :
    public IdealGasMicroThermoTestBase<Scalar,Antioch::NASA7CurveFit<Scalar> >,
    public NASA7ThermoTestHelper<Scalar>
  {
  protected:

    virtual void init_nasa_coeffs()
    {
      // We only need to populate N2, O2, and NO for these tests
      // Populate "lower temp interval first
      {
        // N2 coeffs for 300K-1000K
        std::vector<Scalar> & N2_coeffs = this->_nasa_coeffs_lower_interval["N2"];
        N2_coeffs.resize(7);
        N2_coeffs[0] =  3.298677000E+00L;
        N2_coeffs[1] =  1.408240400E-03L;
        N2_coeffs[2] = -3.963222000E-06L;
        N2_coeffs[3] =  5.641515000E-09L;
        N2_coeffs[4] = -2.444854000E-12L;
        N2_coeffs[5] = -1.020899900E+03L;
        N2_coeffs[6] =  3.950372000E+00L;

        // O2 coeffs for 200K-1000K
        std::vector<Scalar> & O2_coeffs = this->_nasa_coeffs_lower_interval["O2"];
        O2_coeffs.resize(7);
        O2_coeffs[0] =  3.782456360E+00L;
        O2_coeffs[1] = -2.996734160E-03L;
        O2_coeffs[2] =  9.847302010E-06L;
        O2_coeffs[3] = -9.681295090E-09L;
        O2_coeffs[4] =  3.243728370E-12L;
        O2_coeffs[5] = -1.063943560E+03L;
        O2_coeffs[6] =  3.657675730E+00L;

        // NO coeffs for 200K-1000K
        std::vector<Scalar> & NO_coeffs = this->_nasa_coeffs_lower_interval["NO"];
        NO_coeffs.resize(7);
        NO_coeffs[0] =  4.218476300E+00L;
        NO_coeffs[1] = -4.638976000E-03L;
        NO_coeffs[2] =  1.104102200E-05L;
        NO_coeffs[3] = -9.336135400E-09L;
        NO_coeffs[4] =  2.803577000E-12L;
        NO_coeffs[5] =  9.844623000E+03L;
        NO_coeffs[6] =  2.280846400E+00L;
      }

      // Now populate "upper" temp interval first
      {
        // N2 coeffs for 1000K-5000K
        std::vector<Scalar> & N2_coeffs = this->_nasa_coeffs_upper_interval["N2"];
        N2_coeffs.resize(7);
        N2_coeffs[0] =  2.926640000E+00L;
        N2_coeffs[1] =  1.487976800E-03L;
        N2_coeffs[2] = -5.684760000E-07L;
        N2_coeffs[3] =  1.009703800E-10L;
        N2_coeffs[4] = -6.753351000E-15L;
        N2_coeffs[5] = -9.227977000E+02L;
        N2_coeffs[6] =  5.980528000E+00L;

        // O2 coeffs for 1000K-3500K
        std::vector<Scalar> & O2_coeffs = this->_nasa_coeffs_upper_interval["O2"];
        O2_coeffs.resize(7);
        O2_coeffs[0] =  3.282537840E+00L;
        O2_coeffs[1] =  1.483087540E-03L;
        O2_coeffs[2] = -7.579666690E-07L;
        O2_coeffs[3] =  2.094705550E-10L;
        O2_coeffs[4] = -2.167177940E-14L;
        O2_coeffs[5] = -1.088457720E+03L;
        O2_coeffs[6] =  5.453231290E+00L;

        // NO coeffs for 1000K-6000K
        std::vector<Scalar> & NO_coeffs = this->_nasa_coeffs_upper_interval["NO"];
        NO_coeffs.resize(7);
        NO_coeffs[0] =  3.260605600E+00L;
        NO_coeffs[1] =  1.191104300E-03L;
        NO_coeffs[2] = -4.291704800E-07L;
        NO_coeffs[3] =  6.945766900E-11L;
        NO_coeffs[4] = -4.033609900E-15L;
        NO_coeffs[5] =  9.920974600E+03L;
        NO_coeffs[6] =  6.369302700E+00L;
      }

    }

    virtual void parse_nasa_coeffs( Antioch::NASAThermoMixture<Scalar,Antioch::NASA7CurveFit<Scalar> > & nasa_mixture )
    {
      std::string thermo_filename =
        std::string(ANTIOCH_SHARE_XML_INPUT_FILES_SOURCE_PATH)+"gri30.xml";

       Antioch::read_nasa_mixture_data( nasa_mixture, thermo_filename, Antioch::XML );
    }

  };

  template<typename Scalar>
  class IdealGasMicroThermoTestNASA9Base :
    public IdealGasMicroThermoTestBase<Scalar,Antioch::NASA9CurveFit<Scalar> >,
    public NASA9ThermoTestHelper<Scalar>
  {
  protected:

    virtual void init_nasa_coeffs()
    {
      // We only need to populate N2, O2, and NO for these tests
      // Populate "lower temp interval first
      {
        // N2 coeffs for 200K-1000K
        std::vector<Scalar> & N2_coeffs = this->_nasa_coeffs_lower_interval["N2"];
        N2_coeffs.resize(9);
        N2_coeffs[0] =  2.21037122e+04L;
        N2_coeffs[1] = -3.81846145e+02L;
        N2_coeffs[2] =  6.08273815e+00L;
        N2_coeffs[3] = -8.53091381e-03L;
        N2_coeffs[4] =  1.38464610e-05L;
        N2_coeffs[5] = -9.62579293e-09L;
        N2_coeffs[6] =  2.51970560e-12L;
        N2_coeffs[7] =  7.10845911e+02L;
        N2_coeffs[8] = -1.07600320e+01L;

        // O2 coeffs for 200K-1000K
        std::vector<Scalar> & O2_coeffs = this->_nasa_coeffs_lower_interval["O2"];
        O2_coeffs.resize(9);
        O2_coeffs[0] = -3.42556269e+04L;
        O2_coeffs[1] =  4.84699986e+02L;
        O2_coeffs[2] =  1.11901159e+00L;
        O2_coeffs[3] =  4.29388743e-03L;
        O2_coeffs[4] = -6.83627313e-07L;
        O2_coeffs[5] = -2.02337478e-09L;
        O2_coeffs[6] =  1.03904064e-12L;
        O2_coeffs[7] = -3.39145434e+03L;
        O2_coeffs[8] =  1.84969912e+01L;

        // NO coeffs for 200K-1000K
        std::vector<Scalar> & NO_coeffs = this->_nasa_coeffs_lower_interval["NO"];
        NO_coeffs.resize(9);
        NO_coeffs[0] = -1.14391658e+04L;
        NO_coeffs[1] =  1.53646774e+02L;
        NO_coeffs[2] =  3.43146865e+00L;
        NO_coeffs[3] = -2.66859213e-03L;
        NO_coeffs[4] =  8.48139877e-06L;
        NO_coeffs[5] = -7.68511079e-09L;
        NO_coeffs[6] =  2.38679758e-12L;
        NO_coeffs[7] =  9.09794974e+03L;
        NO_coeffs[8] =  6.72872795e+00L;
      }

      // Now populate "upper" temp interval first
      {
        // N2 coeffs for 1000K-6000K
        std::vector<Scalar> & N2_coeffs = this->_nasa_coeffs_upper_interval["N2"];
        N2_coeffs.resize(9);
        N2_coeffs[0] =  5.87709908e+05L;
        N2_coeffs[1] = -2.23924255e+03L;
        N2_coeffs[2] =  6.06694267e+00L;
        N2_coeffs[3] = -6.13965296e-04L;
        N2_coeffs[4] =  1.49179819e-07L;
        N2_coeffs[5] = -1.92309442e-11L;
        N2_coeffs[6] =  1.06194871e-15L;
        N2_coeffs[7] =  1.28320618e+04L;
        N2_coeffs[8] = -1.58663484e+01L;

        // O2 coeffs for 1000K-6000K
        std::vector<Scalar> & O2_coeffs = this->_nasa_coeffs_upper_interval["O2"];
        O2_coeffs.resize(9);
        O2_coeffs[0] = -1.03793994e+06L;
        O2_coeffs[1] =  2.34483275e+03L;
        O2_coeffs[2] =  1.81972949e+00L;
        O2_coeffs[3] =  1.26784887e-03L;
        O2_coeffs[4] = -2.18807142e-07L;
        O2_coeffs[5] =  2.05372411e-11L;
        O2_coeffs[6] = -8.19349062e-16L;
        O2_coeffs[7] = -1.68901253e+04L;
        O2_coeffs[8] =  1.73871835e+01L;

        // NO coeffs for 1000K-6000K
        std::vector<Scalar> & NO_coeffs = this->_nasa_coeffs_upper_interval["NO"];
        NO_coeffs.resize(9);
        NO_coeffs[0] =  2.23903708e+05L;
        NO_coeffs[1] = -1.28965624e+03L;
        NO_coeffs[2] =  5.43394039e+00L;
        NO_coeffs[3] = -3.65605546e-04L;
        NO_coeffs[4] =  9.88101763e-08L;
        NO_coeffs[5] = -1.41608327e-11L;
        NO_coeffs[6] =  9.38021642e-16L;
        NO_coeffs[7] =  1.75029422e+04L;
        NO_coeffs[8] = -8.50169908e+00L;
      }
    }

    virtual void parse_nasa_coeffs( Antioch::NASAThermoMixture<Scalar,Antioch::NASA9CurveFit<Scalar> > & nasa_mixture )
    {
      std::string thermo_filename =
        std::string(ANTIOCH_SHARE_XML_INPUT_FILES_SOURCE_PATH)+"nasa9_thermo.xml";

      Antioch::read_nasa_mixture_data( nasa_mixture, thermo_filename, Antioch::XML );
    }

  };

#define DEFINE_IDEALGASMICROTHERMO_SCALAR_TEST(Classname,BaseClass,Scalar) \
  class Classname : public BaseClass<Scalar>                            \
  {                                                                     \
  public:                                                               \
    CPPUNIT_TEST_SUITE( Classname );                                    \
    CPPUNIT_TEST(test_cv_trans);                                        \
    CPPUNIT_TEST(test_cv_trans_over_R);                                 \
    CPPUNIT_TEST(test_cv_tr);                                           \
    CPPUNIT_TEST(test_cv_tr_over_R);                                    \
    CPPUNIT_TEST(test_cv_rot);                                          \
    CPPUNIT_TEST(test_cv_rot_over_R);                                   \
    CPPUNIT_TEST(test_cv_tr_mix);                                       \
    CPPUNIT_TEST_SUITE_END();                                           \
  }

  DEFINE_IDEALGASMICROTHERMO_SCALAR_TEST(IdealGasMicroThermoNASA7FloatTest,
                                         IdealGasMicroThermoTestNASA7Base,
                                         float);
  DEFINE_IDEALGASMICROTHERMO_SCALAR_TEST(IdealGasMicroThermoNASA7DoubleTest,
                                         IdealGasMicroThermoTestNASA7Base,
                                         double);
  DEFINE_IDEALGASMICROTHERMO_SCALAR_TEST(IdealGasMicroThermoNASA7LongDoubleTest,
                                         IdealGasMicroThermoTestNASA7Base,
                                         long double);

  DEFINE_IDEALGASMICROTHERMO_SCALAR_TEST(IdealGasMicroThermoNASA9FloatTest,
                                         IdealGasMicroThermoTestNASA9Base,
                                         float);
  DEFINE_IDEALGASMICROTHERMO_SCALAR_TEST(IdealGasMicroThermoNASA9DoubleTest,
                                         IdealGasMicroThermoTestNASA9Base,
                                         double);
  DEFINE_IDEALGASMICROTHERMO_SCALAR_TEST(IdealGasMicroThermoNASA9LongDoubleTest,
                                         IdealGasMicroThermoTestNASA9Base,
                                         long double);

  CPPUNIT_TEST_SUITE_REGISTRATION( IdealGasMicroThermoNASA7FloatTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( IdealGasMicroThermoNASA7DoubleTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( IdealGasMicroThermoNASA7LongDoubleTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( IdealGasMicroThermoNASA9FloatTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( IdealGasMicroThermoNASA9DoubleTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( IdealGasMicroThermoNASA9LongDoubleTest );

} // end namespace AntiochTesting

#endif // ANTIOCH_HAVE_CPPUNIT
