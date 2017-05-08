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

// Antioch
#include "antioch/chemical_mixture.h"
#include "antioch/stat_mech_thermo.h"
#include "antioch/ideal_gas_micro_thermo.h"
#include "antioch/nasa_mixture.h"
#include "antioch/nasa_evaluator.h"
#include "antioch/nasa7_curve_fit.h"
#include "antioch/nasa9_curve_fit.h"

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

    void test_cv_trans_over_R()
    {
      for( unsigned int s = 0; s < this->_n_species; s++ )
        this->test_scalar_rel( Scalar(1.5),
                               _thermo->cv_trans_over_R(s),
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
    CPPUNIT_TEST(test_cv_trans_over_R);                         \
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

    void test_cv_trans_over_R()
    {
      for( unsigned int s = 0; s < this->_n_species; s++ )
        this->test_scalar_rel( Scalar(1.5),
                               _thermo->cv_trans_over_R(s),
                               std::numeric_limits<Scalar>::epsilon()*10 );
    }

    void setUp()
    {
      this->init();
      _nasa_mixture = new Antioch::NASAThermoMixture<Scalar,NASACurveFit>( *(this->_chem_mixture) );
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

  private:

    Antioch::NASAThermoMixture<Scalar,NASACurveFit> * _nasa_mixture;

    Antioch::NASAEvaluator<Scalar,NASACurveFit> * _nasa_evaluator;

    Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<Scalar,NASACurveFit>,Scalar> * _thermo;

  };

#define DEFINE_IDEALGASMICROTHERMO_SCALAR_TEST(classname,scalar,CurveFit) \
  class classname :                      \
    public IdealGasMicroThermoTestBase<scalar,Antioch::CurveFit<scalar> > \
  {                                                                     \
  public:                                                               \
    CPPUNIT_TEST_SUITE( classname );     \
    CPPUNIT_TEST(test_cv_trans_over_R);                                 \
    CPPUNIT_TEST_SUITE_END();                                           \
  }

  DEFINE_IDEALGASMICROTHERMO_SCALAR_TEST(IdealGasMicroThermoNASA7FloatTest,float,NASA7CurveFit);
  DEFINE_IDEALGASMICROTHERMO_SCALAR_TEST(IdealGasMicroThermoNASA7DoubleTest,double,NASA7CurveFit);
  DEFINE_IDEALGASMICROTHERMO_SCALAR_TEST(IdealGasMicroThermoNASA7LongDoubleTest,long double,NASA7CurveFit);
  DEFINE_IDEALGASMICROTHERMO_SCALAR_TEST(IdealGasMicroThermoNASA9FloatTest,float,NASA9CurveFit);
  DEFINE_IDEALGASMICROTHERMO_SCALAR_TEST(IdealGasMicroThermoNASA9DoubleTest,double,NASA9CurveFit);
  DEFINE_IDEALGASMICROTHERMO_SCALAR_TEST(IdealGasMicroThermoNASA9LongDoubleTest,long double,NASA9CurveFit);

  CPPUNIT_TEST_SUITE_REGISTRATION( IdealGasMicroThermoNASA7FloatTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( IdealGasMicroThermoNASA7DoubleTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( IdealGasMicroThermoNASA7LongDoubleTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( IdealGasMicroThermoNASA9FloatTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( IdealGasMicroThermoNASA9DoubleTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( IdealGasMicroThermoNASA9LongDoubleTest );

} // end namespace AntiochTesting

#endif // ANTIOCH_HAVE_CPPUNIT
