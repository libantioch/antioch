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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "antioch_config.h"

#include <valarray>

#ifdef ANTIOCH_HAVE_EIGEN
#include "Eigen/Dense"
#endif

#ifdef ANTIOCH_HAVE_METAPHYSICL
#include "metaphysicl/numberarray.h"
#endif

#ifdef ANTIOCH_HAVE_VEXCL
#include "vexcl/vexcl.hpp"
#endif

// Antioch
// Declare metaprogramming overloads before they're used
#include "antioch/eigen_utils_decl.h"
#include "antioch/metaphysicl_utils_decl.h"
#include "antioch/valarray_utils_decl.h"
#include "antioch/vexcl_utils_decl.h"

#include "antioch/chemical_mixture.h"
#include "antioch/cea_evaluator.h"
#include "antioch/physical_constants.h"
#include "antioch/cea_mixture_ascii_parsing.h"

#include "antioch/eigen_utils.h"
#include "antioch/metaphysicl_utils.h"
#include "antioch/valarray_utils.h"
#include "antioch/vexcl_utils.h"

#ifdef ANTIOCH_HAVE_GRVY
#include "grvy.h"

GRVY::GRVY_Timer_Class gt;
#endif

// C++
#include <cmath>
#include <limits>

template <typename Scalar, typename TrioScalars>
int test_cp( const std::string& species_name, unsigned int species,
	     TrioScalars cp_exact, TrioScalars T,
	     const Antioch::CEAEvaluator<Scalar>& thermo,
             const std::string& testname )
{
  using std::abs;

  int return_flag = 0;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 25;

  typedef typename Antioch::
		     template TempCache<TrioScalars> Cache;

#ifdef ANTIOCH_HAVE_GRVY
  gt.BeginTimer(testname);
#endif

  const TrioScalars cp = thermo.cp(Cache(T), species);

#ifdef ANTIOCH_HAVE_GRVY
  gt.EndTimer(testname);
#endif

  // Workaround for a non-standard gcc definition of valarray binary operator-
  const TrioScalars diff = cp_exact - cp; 
  const TrioScalars rel_cp_error = abs(diff/cp_exact);

  if( Antioch::max(rel_cp_error) > tol )
    {
      std::cerr << "Error: Mismatch in species specific heat."
		<< std::setprecision
		     (std::numeric_limits<Scalar>::digits10 + 1)
		<< "\nspecies    = " << species_name
		<< "\ncp         = " << cp
		<< "\ncp_exact   = " << cp_exact
		<< "\ndifference = " << diff
		<< "\nrelative   = " << rel_cp_error
		<< "\ntolerance  = " << tol
		<< "\nT = " << T << std::endl;
      return_flag = 1;
    }

  return return_flag;
}


template <typename Scalar>
Scalar cp( Scalar T, Scalar a0, Scalar a1, Scalar a2, 
	   Scalar a3, Scalar a4, Scalar a5, Scalar a6 )
{
  if( T < 200.1)
    T = 200.1;

  return a0/(T*T) + a1/T + a2 + a3*T + a4*(T*T) + a5*(T*T*T) + a6*(T*T*T*T);
}


template <typename TrioScalars>
int vectester(const TrioScalars& example, const std::string& testname)
{
  typedef typename Antioch::value_type<TrioScalars>::type Scalar;

  std::vector<std::string> species_str_list;
  const unsigned int n_species = 5;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "O2" );
  species_str_list.push_back( "N" );
  species_str_list.push_back( "O" );
  species_str_list.push_back( "NO" );

  const Scalar Mm_N = 14.008e-3L;
  const Scalar Mm_O = 16.000e-3L;
  const Scalar Mm_N2 = 2.L * Mm_N;
  const Scalar Mm_O2 = 2.L * Mm_O;
  const Scalar Mm_NO = Mm_N + Mm_O;

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );

  Antioch::CEAThermoMixture<Scalar> cea_mixture( chem_mixture );
  Antioch::read_cea_mixture_data_ascii_default( cea_mixture );
  Antioch::CEAEvaluator<Scalar> thermo( cea_mixture );

  //const Scalar P = 100000.0;
  TrioScalars T = example;
  for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
    {
      T[3*tuple] = 190.0;
      T[3*tuple+1] = 1500.0;
      T[3*tuple+2] = 10000.0;
    }

  const Scalar R_N2 = Antioch::Constants::R_universal<Scalar>()/Mm_N2;
  const Scalar R_O2 = Antioch::Constants::R_universal<Scalar>()/Mm_O2;
  const Scalar R_N = Antioch::Constants::R_universal<Scalar>()/Mm_N;
  const Scalar R_O = Antioch::Constants::R_universal<Scalar>()/Mm_O;
  const Scalar R_NO = Antioch::Constants::R_universal<Scalar>()/Mm_NO;

  int return_flag = 0;

  // Test N2 cp
  {
    unsigned int index = 0;
    TrioScalars cp_N2 = example;
    for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
      {
        cp_N2[3*tuple  ] = R_N2*cp( Scalar(T[0]), Scalar(2.21037122e+04), Scalar(-3.81846145e+02), Scalar(6.08273815e+00), 
			    Scalar(-8.53091381e-03),  Scalar(1.38464610e-05), Scalar(-9.62579293e-09),  Scalar(2.51970560e-12));

        cp_N2[3*tuple+1] = R_N2*cp( Scalar(T[1]), Scalar(5.87709908e+05), Scalar(-2.23924255e+03),  Scalar(6.06694267e+00),
			    Scalar(-6.13965296e-04), Scalar(1.49179819e-07), Scalar(-1.92309442e-11), Scalar(1.06194871e-15) );

        cp_N2[3*tuple+2] = R_N2*cp( Scalar(T[2]), Scalar(8.30971200e+08), Scalar(-6.42048187e+05),  Scalar(2.02020507e+02), Scalar(-3.06501961e-02),
			Scalar(2.48685558e-06), Scalar(-9.70579208e-11), Scalar(1.43751673e-15));
      }
 
    const Antioch::Species species = chem_mixture.species_list()[index];
    const std::string species_name = chem_mixture.species_inverse_name_map().find(species)->second;

    int return_flag_temp = 0;
    return_flag_temp = test_cp( species_name, index, cp_N2, T, thermo, testname );
    if( return_flag_temp != 0 ) return_flag = 1;
  }

  // Test O2 cp
  {
    unsigned int index = 1;
    TrioScalars cp_O2 = example;
    for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
      {
        cp_O2[3*tuple  ] = R_O2*cp( Scalar(T[0]), Scalar(-3.42556269e+04), Scalar(4.84699986e+02), Scalar(1.11901159e+00), 
			    Scalar(4.29388743e-03), Scalar(-6.83627313e-07), Scalar(-2.02337478e-09) , Scalar(1.03904064e-12) );

        cp_O2[3*tuple+1] = R_O2*cp( Scalar(T[1]), Scalar(-1.03793994e+06), Scalar(2.34483275e+03), Scalar(1.81972949e+00), Scalar(1.26784887e-03), 
			    Scalar(-2.18807142e-07), Scalar(2.05372411e-11), Scalar(-8.19349062e-16) );

        cp_O2[3*tuple+2] = R_O2*cp( Scalar(T[2]), Scalar(4.97515261e+08), Scalar(-2.86602339e+05), Scalar(6.69015464e+01), Scalar(-6.16971869e-03),  
			    Scalar(3.01623757e-07), Scalar(-7.42087888e-12), Scalar(7.27744063e-17));
      }
 
    const Antioch::Species species = chem_mixture.species_list()[index];
    const std::string species_name = chem_mixture.species_inverse_name_map().find(species)->second;

    int return_flag_temp = 0;
    return_flag_temp = test_cp( species_name, index, cp_O2, T, thermo, testname );
    if( return_flag_temp != 0 ) return_flag = 1;
  }

  // Test N cp
  {
    unsigned int index = 2;
    TrioScalars cp_N = example;
    for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
      {
        cp_N[3*tuple  ] = R_N*cp( Scalar(T[0]), Scalar(0.00000000e+00), Scalar(0.00000000e+00), Scalar(2.50000000e+00), Scalar(0.00000000e+00),
		          Scalar(0.00000000e+00), Scalar(0.00000000e+00), Scalar(0.00000000e+00));

        cp_N[3*tuple+1] = R_N*cp( Scalar(T[1]), Scalar(8.87650138e+04), Scalar(-1.07123150e+02), Scalar(2.36218829e+00), Scalar(2.91672008e-04),
		          Scalar(-1.72951510e-07), Scalar(4.01265788e-11), Scalar(-2.67722757e-15) );

        cp_N[3*tuple+2] = R_N*cp( Scalar(T[2]), Scalar(5.47518105e+08), Scalar(-3.10757498e+05), Scalar(6.91678274e+01), Scalar(-6.84798813e-03),
		          Scalar(3.82757240e-07), Scalar(-1.09836771e-11), Scalar(1.27798602e-16));
      }
 
    const Antioch::Species species = chem_mixture.species_list()[index];
    const std::string species_name = chem_mixture.species_inverse_name_map().find(species)->second;

    int return_flag_temp = 0;
    return_flag_temp = test_cp( species_name, index, cp_N, T, thermo, testname );
    if( return_flag_temp != 0 ) return_flag = 1;
  }


  // Test O cp
  {
    unsigned int index = 3;
    TrioScalars cp_O = example;
    for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
      {
        cp_O[3*tuple  ] = R_O*cp( Scalar(T[0]), Scalar(-7.95361130e+03) , Scalar(1.60717779e+02), Scalar(1.96622644e+00), Scalar(1.01367031e-03),
		          Scalar(-1.11041542e-06), Scalar(6.51750750e-10), Scalar(-1.58477925e-13) );

        cp_O[3*tuple+1] = R_O*cp( Scalar(T[1]), Scalar(2.61902026e+05), Scalar(-7.29872203e+02), Scalar(3.31717727e+00), Scalar(-4.28133436e-04), 
		          Scalar(1.03610459e-07), Scalar(-9.43830433e-12), Scalar(2.72503830e-16) );

        cp_O[3*tuple+2] = R_O*cp( Scalar(T[2]), Scalar(1.77900426e+08), Scalar(-1.08232826e+05), Scalar(2.81077837e+01), Scalar(-2.97523226e-03),
		          Scalar(1.85499753e-07), Scalar(-5.79623154e-12), Scalar(7.19172016e-17) );
      }
 
    const Antioch::Species species = chem_mixture.species_list()[index];
    const std::string species_name = chem_mixture.species_inverse_name_map().find(species)->second;

    int return_flag_temp = 0;
    return_flag_temp = test_cp( species_name, index, cp_O, T, thermo, testname );
    if( return_flag_temp != 0 ) return_flag = 1;
  }


  // Test NO cp
  {
    unsigned int index = 4;
    TrioScalars cp_NO = example;
    for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
      {
        cp_NO[3*tuple  ] = R_NO*cp( Scalar(T[0]), Scalar(-1.14391658e+04), Scalar(1.53646774e+02), Scalar(3.43146865e+00), Scalar(-2.66859213e-03),
			    Scalar(8.48139877e-06), Scalar(-7.68511079e-09), Scalar(2.38679758e-12) );

        cp_NO[3*tuple+1] = R_NO*cp( Scalar(T[1]), Scalar(2.23903708e+05), Scalar(-1.28965624e+03), Scalar(5.43394039e+00), Scalar(-3.65605546e-04), 
			    Scalar(9.88101763e-08), Scalar(-1.41608327e-11), Scalar(9.38021642e-16) );

        cp_NO[3*tuple+2] = R_NO*cp( Scalar(T[2]), Scalar(-9.57530764e+08), Scalar(5.91243671e+05), Scalar(-1.38456733e+02), Scalar(1.69433998e-02), 
			    Scalar(-1.00735146e-06), Scalar(2.91258526e-11), Scalar(-3.29511091e-16) );
      }

    const Antioch::Species species = chem_mixture.species_list()[index];
    const std::string species_name = chem_mixture.species_inverse_name_map().find(species)->second;

    int return_flag_temp = 0;
    return_flag_temp = test_cp( species_name, index, cp_NO, T, thermo, testname );
    if( return_flag_temp != 0 ) return_flag = 1;
  }

  return return_flag;
}


int main()
{
  int returnval = 0;

  returnval = returnval ||
    vectester (std::valarray<float>(3*ANTIOCH_N_TUPLES), "valarray<float>");
  returnval = returnval ||
    vectester (std::valarray<double>(3*ANTIOCH_N_TUPLES), "valarray<double>");
// We're not getting the full long double precision yet?
//  returnval = returnval ||
//    vectester<long double, std::valarray<long double> >
//      (std::valarray<long double>(3*ANTIOCH_N_TUPLES), "valarray<ld>");
#ifdef ANTIOCH_HAVE_EIGEN
  returnval = returnval ||
    vectester (Eigen::Array<float, 3*ANTIOCH_N_TUPLES, 1>(), "Eigen::ArrayXf");
  returnval = returnval ||
    vectester (Eigen::Array<double, 3*ANTIOCH_N_TUPLES, 1>(), "Eigen::ArrayXd");
//  returnval = returnval ||
//    vectester (Eigen::Array<long double, 3*ANTIOCH_N_TUPLES, 1>(), "Eigen::ArrayXld");
#endif
#ifdef ANTIOCH_HAVE_METAPHYSICL
  returnval = returnval ||
    vectester (MetaPhysicL::NumberArray<3*ANTIOCH_N_TUPLES, float> (0), "NumberArray<float>");
  returnval = returnval ||
    vectester (MetaPhysicL::NumberArray<3*ANTIOCH_N_TUPLES, double> (0), "NumberArray<double>");
//  returnval = returnval ||
//    vectester (MetaPhysicL::NumberArray<3*ANTIOCH_N_TUPLES, long double> (0)), "NumberArray<ld>");
#endif
#ifdef ANTIOCH_HAVE_VEXCL
  vex::Context ctx_f (vex::Filter::All);
  if (!ctx_f.empty())
    returnval = returnval ||
      vectester (vex::vector<float> (ctx_f, 3*ANTIOCH_N_TUPLES), "vex::vector<float>");

  vex::Context ctx_d (vex::Filter::DoublePrecision);
  if (!ctx_d.empty())
    returnval = returnval ||
      vectester (vex::vector<double> (ctx_d, 3*ANTIOCH_N_TUPLES), "vex::vector<double>");
#endif

#ifdef ANTIOCH_HAVE_GRVY
  gt.Finalize();
  gt.Summarize();
#endif

  return returnval;
}
