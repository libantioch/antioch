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

#ifndef ANTIOCH_NASA7_MIXTURE_VECTOR_TEST_BASE_H
#define ANTIOCH_NASA7_MIXTURE_VECTOR_TEST_BASE_H

#include "antioch_config.h"

#ifdef ANTIOCH_HAVE_CPPUNIT

#include "thermo_vector_test_base.h"
#include "nasa7_thermo_test_base.h"

// Antioch
#include "antioch/chemical_mixture.h"
#include "antioch/nasa7_curve_fit.h"
#include "antioch/nasa_mixture.h"
#include "antioch/nasa_mixture_parsing.h"
#include "antioch/xml_parser.h"

namespace AntiochTesting
{
  template<typename PairScalars>
  class NASA7MixtureVectorTestBase :
    public ThermoVectorTestBase<Antioch::NASA7CurveFit<typename Antioch::value_type<PairScalars>::type>, PairScalars>,
    public NASA7ThermoTestBase<typename Antioch::value_type<PairScalars>::type>
  {

  private:

    typedef typename Antioch::value_type<PairScalars>::type Scalar;
    typedef typename Antioch::NASA7CurveFit<Scalar> NASA7Fit;

  public:

    virtual void init()
    {
      NASA7ThermoTestBase<Scalar>::init();

      std::vector<std::string> species_str_list(2);
      species_str_list[0] = "H2";
      species_str_list[1] = "N2";

      this->_chem_mixture = new Antioch::ChemicalMixture<Scalar> ( species_str_list );

      std::string thermo_filename = std::string(ANTIOCH_TESTING_INPUT_FILES_PATH)+"nasa7_thermo_test.xml";

      this->_nasa_mixture = new Antioch::NASAThermoMixture<Scalar,NASA7Fit>( *(this->_chem_mixture) );

      Antioch::read_nasa_mixture_data( *(this->_nasa_mixture), thermo_filename, Antioch::XML );
    }

    virtual void clear()
    {
      delete this->_nasa_mixture;
      delete this->_chem_mixture;
    }

  protected:

    virtual void prep_coeffs( Scalar T, unsigned int species,
                              std::vector<Scalar>& coeffs )
    {
      unsigned int n_coeffs = 7;
      coeffs.resize(n_coeffs);

      std::string species_name =
        this->_chem_mixture->chemical_species()[species]->species();

      if( species == 0 )
        {
          CPPUNIT_ASSERT_EQUAL( std::string("H2"), species_name );

          if( T > 200 && T <= 1000 )
            {
              for( unsigned int i = 0; i < n_coeffs; i++ )
                coeffs[i] = this->_H2_coeffs_200_1000[i];
            }
          else if( T > 1000 && T <= 3500 )
            {
              for( unsigned int i = 0; i < n_coeffs; i++ )
                coeffs[i] = this->_H2_coeffs_1000_3500[i];
            }
          else
            CPPUNIT_FAIL("ERROR: Invalid temperature for species H2!");
        }
      else if( species == 1 )
        {
          CPPUNIT_ASSERT_EQUAL( std::string("N2"), species_name );

          if( T > 300 && T <= 1000 )
            {
              for( unsigned int i = 0; i < n_coeffs; i++ )
                coeffs[i] = this->_N2_coeffs_300_1000[i];
            }
          else if( T > 1000 && T <= 5000 )
            {
              for( unsigned int i = 0; i < n_coeffs; i++ )
                coeffs[i] = this->_N2_coeffs_1000_5000[i];
            }
          else
            CPPUNIT_FAIL("ERROR: Invalid temperature for species N2!");
        }
      else
        CPPUNIT_FAIL("ERROR: Invalid species index!");
    }

    virtual Scalar cp_over_R_exact( Scalar T, const std::vector<Scalar>& coeffs )
    {
      CPPUNIT_ASSERT_MESSAGE("ERROR: Not enough coeffs!",
                             coeffs.size() >= 5 );
      return this->cp_exact( T,
                             coeffs[0],
                             coeffs[1],
                             coeffs[2],
                             coeffs[3],
                             coeffs[4] );
    }

    virtual Scalar h_over_RT_exact( Scalar T, const std::vector<Scalar>& coeffs )
    {
      CPPUNIT_ASSERT_MESSAGE("ERROR: Not enough coeffs!",
                             coeffs.size() >= 6 );

      return this->h_exact( T,
                            coeffs[0],
                            coeffs[1],
                            coeffs[2],
                            coeffs[3],
                            coeffs[4],
                            coeffs[5] );
    }

    virtual Scalar s_over_R_exact( Scalar T, const std::vector<Scalar>& coeffs )
    {
      CPPUNIT_ASSERT_MESSAGE("ERROR: Not enough coeffs!",
                             coeffs.size() >= 7 );

      return this->s_exact( T,
                            coeffs[0],
                            coeffs[1],
                            coeffs[2],
                            coeffs[3],
                            coeffs[4],
                            coeffs[6] );
    }

  };

} // end namespace AntiochTesting

#endif // ANTIOCH_HAVE_CPPUNIT

#endif // ANTIOCH_NASA7_MIXTURE_VECTOR_TEST_BASE_H
