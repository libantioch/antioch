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

#ifndef ANTIOCH_THERMO_VECTOR_TEST_BASE_H
#define ANTIOCH_THERMO_VECTOR_TEST_BASE_H

#include "antioch_config.h"

#ifdef ANTIOCH_HAVE_CPPUNIT

// CppUnit
#include <cppunit/extensions/HelperMacros.h>

// Antioch
#include "antioch/nasa_evaluator.h"
#include "antioch/temp_cache.h"

namespace AntiochTesting
{
  template <typename NASAFit, typename PairScalars>
  class ThermoVectorTestBase
  {
  private:

    typedef typename Antioch::value_type<PairScalars>::type Scalar;

  public:

    Scalar tol()
    { return std::numeric_limits<Scalar>::epsilon() * 500; }

    //! Test that we actually run from the CppUnit driver
    void test_all_cp()
    {
      Antioch::NASAEvaluator<Scalar,NASAFit> thermo_evaluator( (*_nasa_mixture) );
      PairScalars T = this->setup_T(*(this->_example));

      // Test for each species (0 or 1, species names depend on subclass)
      this->test_cp(thermo_evaluator,0,T,this->tol());
      this->test_cp(thermo_evaluator,1,T,this->tol());
    }

    //! Test that we actually run from the CppUnit driver
    void test_all_h()
    {
      Antioch::NASAEvaluator<Scalar,NASAFit> thermo_evaluator( (*_nasa_mixture) );
      PairScalars T = this->setup_T(*(this->_example));

      // Test for each species (0 or 1, species names depend on subclass)
      this->test_h(thermo_evaluator,0,T,this->tol());
      this->test_h(thermo_evaluator,1,T,this->tol());
    }

    //! Test that we actually run from the CppUnit driver
    void test_all_s()
    {
      Antioch::NASAEvaluator<Scalar,NASAFit> thermo_evaluator( (*_nasa_mixture) );
      PairScalars T = this->setup_T(*(this->_example));

      // Test for each species (0 or 1, species names depend on subclass)
      this->test_s(thermo_evaluator,0,T,this->tol());
      this->test_s(thermo_evaluator,1,T,this->tol());
    }

  protected:

    Antioch::ChemicalMixture<Scalar>* _chem_mixture;
    Antioch::NASAThermoMixture<Scalar,NASAFit>* _nasa_mixture;

    void test_cp( const Antioch::NASAEvaluator<Scalar,NASAFit>& thermo,
                  unsigned int species,
                  const PairScalars& T,
                  Scalar tol )
    {
      Antioch::TempCache<PairScalars> cache(T);
      const PairScalars cp = thermo.cp_over_R(cache,species);

      const PairScalars exact_cp = this->exact_cp_over_R(cache,species);

      for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL( cp[2*tuple],
                                        exact_cp[0],
                                        tol );

          CPPUNIT_ASSERT_DOUBLES_EQUAL( cp[2*tuple+1],
                                        exact_cp[1],
                                        tol );
        }
    }

    void test_h( const Antioch::NASAEvaluator<Scalar,NASAFit>& thermo,
                 unsigned int species,
                 const PairScalars& T,
                 Scalar tol )
    {
      Antioch::TempCache<PairScalars> cache(T);
      const PairScalars h = thermo.h_over_RT(cache,species);

      const PairScalars exact_h = this->exact_h_over_RT(cache,species);

      for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL( h[2*tuple],
                                        exact_h[0],
                                        tol );

          CPPUNIT_ASSERT_DOUBLES_EQUAL( h[2*tuple+1],
                                        exact_h[1],
                                        tol );
        }
    }

    void test_s( const Antioch::NASAEvaluator<Scalar,NASAFit>& thermo,
                 unsigned int species,
                 const PairScalars& T,
                 Scalar tol )
    {
      Antioch::TempCache<PairScalars> cache(T);
      const PairScalars s = thermo.s_over_R(cache,species);

      const PairScalars exact_s = this->exact_s_over_R(cache,species);

      for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL( s[2*tuple],
                                        exact_s[0],
                                        tol );

          CPPUNIT_ASSERT_DOUBLES_EQUAL( s[2*tuple+1],
                                        exact_s[1],
                                        tol );
        }
    }

    // Should be new'd/deleted in subclasses for each PairScalar type
    PairScalars* _example;

    PairScalars setup_T( const PairScalars& example )
    {
      // Construct from example to avoid resizing issues
      PairScalars T = example;
      for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
        {
          // Make sure we hit two different intervals
          T[2*tuple]   = 500.1;
          T[2*tuple+1] = 1600.1;
        }
      return T;
    }

    PairScalars exact_cp_over_R( const Antioch::TempCache<PairScalars>& cache, unsigned species )
    {
      PairScalars e_cp = *(this->_example);
      for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
        {
          e_cp[2*tuple]   =  this->eval_cp_over_R( cache.T[2*tuple], species );
          e_cp[2*tuple+1] =  this->eval_cp_over_R( cache.T[2*tuple+1], species );
        }
      return e_cp;
    }

    PairScalars exact_s_over_R( const Antioch::TempCache<PairScalars>& cache, unsigned species )
    {
      PairScalars e_s = *(this->_example);
      for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
        {
          e_s[2*tuple]   =  this->eval_s_over_R( cache.T[2*tuple], species );
          e_s[2*tuple+1] =  this->eval_s_over_R( cache.T[2*tuple+1], species );
        }
      return e_s;
    }

    PairScalars exact_h_over_RT( const Antioch::TempCache<PairScalars>& cache, unsigned species )
    {
      PairScalars e_h = *(this->_example);
      for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
        {
          e_h[2*tuple]   =  this->eval_h_over_RT( cache.T[2*tuple], species );
          e_h[2*tuple+1] =  this->eval_h_over_RT( cache.T[2*tuple+1], species );
        }
      return e_h;
    }

    Scalar eval_cp_over_R( Scalar T, unsigned int species )
    {
      std::vector<Scalar> coeffs;

      this->prep_coeffs(T,species,coeffs);

      return this->cp_over_R_exact(T,coeffs);
    }

    Scalar eval_h_over_RT( Scalar T, unsigned int species )
    {
      std::vector<Scalar> coeffs;

      this->prep_coeffs(T,species,coeffs);

      return this->h_over_RT_exact(T,coeffs);
    }

    Scalar eval_s_over_R( Scalar T, unsigned int species )
    {
      std::vector<Scalar> coeffs;

      this->prep_coeffs(T,species,coeffs);

      return this->s_over_R_exact(T,coeffs);
    }

    virtual void prep_coeffs( Scalar T, unsigned int species,
                              std::vector<Scalar>& coeffs ) =0;

    virtual Scalar cp_over_R_exact( Scalar T, const std::vector<Scalar>& coeffs ) =0;
    virtual Scalar h_over_RT_exact( Scalar T, const std::vector<Scalar>& coeffs ) =0;
    virtual Scalar s_over_R_exact( Scalar T, const std::vector<Scalar>& coeffs ) =0;

  };

} // end namespace AntiochTesting

#endif // ANTIOCH_HAVE_CPPUNIT

#endif // ANTIOCH_THERMO_VECTOR_TEST_BASE_H
