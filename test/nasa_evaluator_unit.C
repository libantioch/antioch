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

// C++
#include <cmath>
#include <limits>

// Antioch
#include "antioch_config.h"
#include "antioch/physical_constants.h"
#include "antioch/chemical_mixture.h"
#include "antioch/cea_curve_fit.h"
#include "antioch/nasa_curve_fit.h"
#include "antioch/default_filename.h"
#include "antioch/temp_cache.h"
#include "antioch/cea_mixture.h"
#include "antioch/cea_evaluator.h"
#include "antioch/cea_mixture_ascii_parsing.h"
#include "antioch/nasa_mixture_ascii_parsing.h"

template <typename Scalar, typename NASAFit>
int test_cp( const std::string& species_name, unsigned int species, Scalar cp_exact, Scalar T,
	     const Antioch::NASAEvaluator<Scalar,NASAFit>& thermo )
{
  using std::abs;

  int return_flag = 0;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 1000;

  typedef typename Antioch::template TempCache<Scalar> Cache;

  const Scalar cp = thermo.cp(Cache(T), species);

  if( abs( (cp_exact - cp)/cp_exact ) > tol )
    {
      std::cerr << std::scientific << std::setprecision(16)
                << "Error: Mismatch in species specific heat."
		<< "\nspecies    = " << species_name
		<< "\ncp         = " << cp
		<< "\ncp_exact   = " << cp_exact
		<< "\ndifference = " << (cp_exact - cp)
		<< "\ntolerance  = " << tol
		<< "\nT = " << T << std::endl;
      return_flag = 1;
    }

  return return_flag;
}

template <typename Scalar, typename NASAFit>
int test_h( const std::string& species_name, unsigned int species, Scalar h_exact, Scalar T,
            const Antioch::NASAEvaluator<Scalar,NASAFit>& thermo )
{
  using std::abs;

  int return_flag = 0;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 10;

  typedef typename Antioch::template TempCache<Scalar> Cache;

  const Scalar h = thermo.h(Cache(T), species);

  if( abs( (h_exact - h)/h_exact ) > tol )
    {
      std::cerr << std::scientific << std::setprecision(16)
                << "Error: Mismatch in species total enthalpy."
		<< "\nspecies    = " << species_name
		<< "\nh          = " << h
		<< "\nh_exact    = " << h_exact
		<< "\ndifference = " << (h_exact - h)
		<< "\ntolerance  = " << tol
		<< "\nT = " << T << std::endl;
      return_flag = 1;
    }

  return return_flag;
}

template <typename Scalar, typename NASAFit>
int test_s( const std::string& species_name, unsigned int species, Scalar s_exact, Scalar T,
	     const Antioch::NASAEvaluator<Scalar, NASAFit>& thermo )
{
  using std::abs;

  int return_flag = 0;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 1000;

  typedef typename Antioch::template TempCache<Scalar> Cache;

  const Scalar s = thermo.s_over_R(Cache(T), species);

  if( abs( (s_exact - s)/s_exact ) > tol )
    {
      std::cerr << std::scientific << std::setprecision(16)
                << "Error: Mismatch in species entropie."
		<< "\nspecies    = " << species_name
		<< "\ns          = " << s
		<< "\ns_exact    = " << s_exact
		<< "\ndifference = " << (s_exact - s)
		<< "\ntolerance  = " << tol
		<< "\nT = " << T << std::endl;
      return_flag = 1;
    }

  return return_flag;
}

template <typename Scalar, typename NASAFit>
int test_g( const std::string& species_name, unsigned int species, Scalar g_exact, Scalar T,
	     const Antioch::NASAEvaluator<Scalar, NASAFit>& thermo )
{
  using std::abs;

  int return_flag = 0;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 1000;

  typedef typename Antioch::template TempCache<Scalar> Cache;

  const Scalar g =   thermo.h_over_RT(Cache(T), species)
                   - thermo.s_over_R(Cache(T), species) ;

  if( abs( (g_exact - g)/g_exact ) > tol )
    {
      std::cerr << std::scientific << std::setprecision(16)
                << "Error: Mismatch in species free energy (Gibbs energy)."
		<< "\nspecies    = " << species_name
		<< "\ng          = " << g
		<< "\ng_exact    = " << g_exact
		<< "\ndifference = " << (g_exact - g)
		<< "\ntolerance  = " << tol
		<< "\nT = " << T << std::endl;
      return_flag = 1;
    }

  return return_flag;
}

template <typename Scalar>
Scalar cea_cp( Scalar T, Scalar a0, Scalar a1, Scalar a2, 
	   Scalar a3, Scalar a4, Scalar a5, Scalar a6 )
{
  if( T < 200.1)
    T = 200.1L;

  return a0/(T*T) + a1/T + a2 + a3*T + a4*(T*T) + a5*(T*T*T) + a6*(T*T*T*T);
}

template <typename Scalar>
Scalar cea_h( Scalar T, Scalar a0, Scalar a1, Scalar a2, 
          Scalar a3, Scalar a4, Scalar a5, Scalar a6,
          Scalar a8 )
{
  return -a0/(T*T) + a1*std::log(T)/T + a2 + a3*T/2.0L + a4*(T*T)/3.0L + a5*(T*T*T)/4.0L + a6*(T*T*T*T)/5.0L + a8/T;
}
template <typename Scalar>
Scalar cea_s( Scalar T, Scalar a0, Scalar a1, Scalar a2, 
          Scalar a3, Scalar a4, Scalar a5, Scalar a6,
          Scalar a9 )
{
  return -a0/(2.L*T*T) - a1/T + a2*std::log(T) + a3*T + a4*(T*T)/2.0L + a5*(T*T*T)/3.0L + a6*(T*T*T*T)/4.0L + a9;
}

template <typename Scalar>
Scalar nasa_cp( Scalar T, Scalar a0, Scalar a1, Scalar a2, 
	   Scalar a3, Scalar a4)
{
  if( T < 200.1)
    T = 200.1L;

  return a0 + a1*T + a2*T*T + a3*T*T*T + a4*T*T*T*T;
}

template <typename Scalar>
Scalar nasa_h( Scalar T, Scalar a0, Scalar a1, Scalar a2, 
          Scalar a3, Scalar a4, Scalar a5)
{
  return  a0 + a1*T/2.0L + a2*T*T/3.0L + a3*T*T*T/4.0L +
          a4*T*T*T*T/5.0L + a5/T;
}
template <typename Scalar>
Scalar nasa_s( Scalar T, Scalar a0, Scalar a1, Scalar a2, 
          Scalar a3, Scalar a4, Scalar a6)
{
  return a0*std::log(T) + a1*T + a2*T*T/2.0L + a3*T*T*T/3.0L +
         a4*T*T*T*T/4.0L + a6;
}



template <typename Scalar>
int tester(const std::string & nasa_filename)
{
  std::vector<std::string> species_str_list;
  const unsigned int n_species = 5;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "O2" );
  species_str_list.push_back( "OH" );
  species_str_list.push_back( "O" );
  species_str_list.push_back( "H2" );

  const Scalar Mm_H =  1.008e-3L;
  const Scalar Mm_N = 14.008e-3L;
  const Scalar Mm_O = 16.000e-3L;
  const Scalar Mm_N2 = 2.L * Mm_N;
  const Scalar Mm_O2 = 2.L * Mm_O;
  const Scalar Mm_H2 = 2.L * Mm_H;
  const Scalar Mm_OH = Mm_H + Mm_O;

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );

// backward compatibility
  Antioch::NASA9ThermoMixture<Scalar> cea_mixture( chem_mixture );
  Antioch::read_cea_mixture_data_ascii( cea_mixture, Antioch::DefaultFilename::thermo_data() );
  Antioch::NASA9Evaluator<Scalar> thermo( cea_mixture );

// explicit
  Antioch::NASAThermoMixture<Scalar, Antioch::NASACurveFit<Scalar> > nasa_mixture( chem_mixture );
  Antioch::read_nasa_mixture_data_ascii( nasa_mixture, nasa_filename );
  Antioch::NASAEvaluator<Scalar, Antioch::NASACurveFit<Scalar> > nasa_thermo( nasa_mixture );

  //const Scalar P = 100000.0;
  const Scalar T1 = 190.0;
  const Scalar T2 = 1500.0;
  const Scalar T3 = 10000.0;
  const Scalar T4 = 800.0;

  const Scalar R_N2 = Antioch::Constants::R_universal<Scalar>()/Mm_N2;
  const Scalar R_O2 = Antioch::Constants::R_universal<Scalar>()/Mm_O2;
  const Scalar R_H2 = Antioch::Constants::R_universal<Scalar>()/Mm_H2;
  const Scalar R_O  = Antioch::Constants::R_universal<Scalar>()/Mm_O;
  const Scalar R_OH = Antioch::Constants::R_universal<Scalar>()/Mm_OH;

// declare everyone here, easier
  const Scalar cea_N2_a1[10] = { 2.21037122e+04, -3.81846145e+02,  6.08273815e+00, -8.53091381e-03,  1.38464610e-05,
                                -9.62579293e-09,  2.51970560e-12,  0.00000000e+00,  7.10845911e+02, -1.07600320e+01};
  const Scalar cea_N2_a2[10] = { 5.87709908e+05, -2.23924255e+03,  6.06694267e+00, -6.13965296e-04,  1.49179819e-07,
                                -1.92309442e-11,  1.06194871e-15,  0.00000000e+00,  1.28320618e+04, -1.58663484e+01};
  const Scalar cea_N2_a3[10] = { 8.30971200e+08, -6.42048187e+05,  2.02020507e+02, -3.06501961e-02,  2.48685558e-06,
                                -9.70579208e-11,  1.43751673e-15,  0.00000000e+00,  4.93850663e+06, -1.67204791e+03};
  const Scalar cea_O2_a1[10] = {-3.42556269e+04,  4.84699986e+02,  1.11901159e+00,  4.29388743e-03, -6.83627313e-07,
                                -2.02337478e-09,  1.03904064e-12,  0.00000000e+00, -3.39145434e+03,  1.84969912e+01};
  const Scalar cea_O2_a2[10] = {-1.03793994e+06,  2.34483275e+03,  1.81972949e+00,  1.26784887e-03, -2.18807142e-07,
                                 2.05372411e-11, -8.19349062e-16,  0.00000000e+00, -1.68901253e+04,  1.73871835e+01};
  const Scalar cea_O2_a3[10] = { 4.97515261e+08, -2.86602339e+05,  6.69015464e+01, -6.16971869e-03,  3.01623757e-07,
                                -7.42087888e-12,  7.27744063e-17,  0.00000000e+00,  2.29348755e+06, -5.53044968e+02};
  const Scalar cea_OH_a1[10] = {-1.99886366e+03,  9.30014295e+01,  3.05085385e+00,  1.52953035e-03, -3.15789256e-06,
                                 3.31544734e-09, -1.13876303e-12,  0.00000000e+00,  3.24004397e+03,  4.67411290e+00};
  const Scalar cea_OH_a2[10] = { 1.01739349e+06, -2.50995748e+03,  5.11654809e+00,  1.30529875e-04, -8.28431902e-08,
                                 2.00647550e-11, -1.55699342e-15,  0.00000000e+00,  2.04452334e+04, -1.10128250e+01};
  const Scalar cea_OH_a3[10] = { 2.84724229e+08, -1.85953704e+05,  5.00825057e+01, -5.14238586e-03,  2.87554326e-07,
                                -8.22883906e-12,  9.56725603e-17,  0.00000000e+00,  1.46864630e+06, -4.02356407e+02};
  const Scalar cea_O_a1[10]  = {-7.95361130e+03,  1.60717779e+02,  1.96622644e+00,  1.01367031e-03, -1.11041542e-06,
                                 6.51750750e-10, -1.58477925e-13,  0.00000000e+00,  2.84036244e+04,  8.40424182e+00};
  const Scalar cea_O_a2[10]  = { 2.61902026e+05, -7.29872203e+02,  3.31717727e+00, -4.28133436e-04,  1.03610459e-07,
                                -9.43830433e-12,  2.72503830e-16,  0.00000000e+00,  3.39242806e+04, -6.67958535e-01};
  const Scalar cea_O_a3[10]  = { 1.77900426e+08, -1.08232826e+05,  2.81077837e+01, -2.97523226e-03,  1.85499753e-07,
                                -5.79623154e-12,  7.19172016e-17,  0.00000000e+00,  8.89094263e+05, -2.18172815e+02};
  const Scalar cea_H2_a1[10] = { 4.07832281e+04, -8.00918545e+02,  8.21470167e+00, -1.26971436e-02,  1.75360493e-05,
                                -1.20286016e-08,  3.36809316e-12,  0.00000000e+00,  2.68248438e+03, -3.04378866e+01};
  const Scalar cea_H2_a2[10] = { 5.60812338e+05, -8.37149134e+02,  2.97536304e+00,  1.25224993e-03, -3.74071842e-07,
                                 5.93662825e-11, -3.60699573e-15,  0.00000000e+00,  5.33981585e+03, -2.20276405e+00};
  const Scalar cea_H2_a3[10] = { 4.96671613e+08, -3.14744812e+05,  7.98388750e+01, -8.41450419e-03,  4.75306044e-07,
                                -1.37180973e-11,  1.60537460e-16,  0.00000000e+00,  2.48835466e+06, -6.69552419e+02};

  const Scalar nasa_N2_a1[7] = { 0.02926640E+02,  0.14879768E-02, -0.05684760E-05,  0.10097038E-09, -0.06753351E-13,
                                -0.09227977E+04,  0.05980528E+02};
  const Scalar nasa_N2_a2[7] = { 0.03298677E+02,  0.14082404E-02, -0.03963222E-04,  0.05641515E-07, -0.02444854E-10,
                                -0.10208999E+04,  0.03950372E+02};
  const Scalar nasa_O2_a1[7] = { 3.28253784E+00,  1.48308754E-03, -7.57966669E-07,  2.09470555E-10, -2.16717794E-14,
                                -1.08845772E+03,  5.45323129E+00}; 
  const Scalar nasa_O2_a2[7] = { 3.78245636E+00, -2.99673416E-03,  9.84730201E-06, -9.68129509E-09,  3.24372837E-12,
                                -1.06394356E+03,  3.65767573E+00};
  const Scalar nasa_OH_a1[7] = { 3.09288767E+00,  5.48429716E-04,  1.26505228E-07, -8.79461556E-11,  1.17412376E-14,
                                 3.85865700E+03,  4.47669610E+00};
  const Scalar nasa_OH_a2[7] = { 3.99201543E+00, -2.40131752E-03,  4.61793841E-06, -3.88113333E-09,  1.36411470E-12,
                                 3.61508056E+03, -1.03925458E-01};
  const Scalar nasa_O_a1[7]  = { 2.56942078E+00, -8.59741137E-05,  4.19484589E-08, -1.00177799E-11,  1.22833691E-15,
                                 2.92175791E+04,  4.78433864E+00};
  const Scalar nasa_O_a2[7]  = { 3.16826710E+00, -3.27931884E-03,  6.64306396E-06, -6.12806624E-09,  2.11265971E-12,
                                 2.91222592E+04,  2.05193346E+00};
  const Scalar nasa_H2_a1[7] = { 3.33727920E+00, -4.94024731E-05,  4.99456778E-07, -1.79566394E-10,  2.00255376E-14,
                                -9.50158922E+02, -3.20502331E+00};
  const Scalar nasa_H2_a2[7] = { 2.34433112E+00,  7.98052075E-03, -1.94781510E-05,  2.01572094E-08, -7.37611761E-12,
                                -9.17935173E+02,  6.83010238E-01};

  int return_flag = 0;


  // Test N2 cp
  {
    unsigned int index = 0;

// CEA
    const Scalar cea_cp_N2_1 = R_N2*cea_cp( T1, cea_N2_a1[0], cea_N2_a1[1], cea_N2_a1[2], cea_N2_a1[3], cea_N2_a1[4], cea_N2_a1[5], cea_N2_a1[6]);

    const Scalar cea_cp_N2_2 = R_N2*cea_cp( T2, cea_N2_a2[0], cea_N2_a2[1], cea_N2_a2[2], cea_N2_a2[3], cea_N2_a2[4], cea_N2_a2[5], cea_N2_a2[6]);

    const Scalar cea_cp_N2_3 = R_N2*cea_cp( T3, cea_N2_a3[0], cea_N2_a3[1], cea_N2_a3[2], cea_N2_a3[3], cea_N2_a3[4], cea_N2_a3[5], cea_N2_a3[6]);
    
    const Antioch::Species species = chem_mixture.species_list()[index];
    const std::string species_name = chem_mixture.species_inverse_name_map().find(species)->second;

    int return_flag_temp = 0;

    return_flag_temp = test_cp( species_name, index, cea_cp_N2_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_cp( species_name, index, cea_cp_N2_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_cp( species_name, index, cea_cp_N2_3, T3, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;


// NASA [300 - 1000] [1000 - 5000]

    const Scalar nasa_cp_N2_1 = R_N2 * nasa_cp(T4, nasa_N2_a1[0], nasa_N2_a1[1], nasa_N2_a1[2], nasa_N2_a1[3], nasa_N2_a1[4]);

    const Scalar nasa_cp_N2_2 = R_N2 * nasa_cp(T2, nasa_N2_a2[0], nasa_N2_a2[1], nasa_N2_a2[2], nasa_N2_a2[3], nasa_N2_a2[4]);

    return_flag_temp = test_cp( species_name, index, nasa_cp_N2_1, T4, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_cp( species_name, index, nasa_cp_N2_2, T2, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

  }

  // Test O2 cp
  {
    unsigned int index = 1;
    const Scalar cp_1 = R_O2*cea_cp( T1, cea_O2_a1[0], cea_O2_a1[1], cea_O2_a1[2], cea_O2_a1[3], cea_O2_a1[4], cea_O2_a1[5], cea_O2_a1[6]);

    const Scalar cp_2 = R_O2*cea_cp( T2, cea_O2_a2[0], cea_O2_a2[1], cea_O2_a2[2], cea_O2_a2[3], cea_O2_a2[4], cea_O2_a2[5], cea_O2_a2[6]);

    const Scalar cp_3 = R_O2*cea_cp( T3, cea_O2_a3[0], cea_O2_a3[1], cea_O2_a3[2], cea_O2_a3[3], cea_O2_a3[4], cea_O2_a3[5], cea_O2_a3[6]);
    
    const Antioch::Species species = chem_mixture.species_list()[index];
    const std::string species_name = chem_mixture.species_inverse_name_map().find(species)->second;

    int return_flag_temp = 0;
    return_flag_temp = test_cp( species_name, index, cp_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_cp( species_name, index, cp_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_cp( species_name, index, cp_3, T3, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

// NASA [300 - 1000] [1000 - 5000]
    const Scalar nasa_cp_O2_1 = R_O2 * nasa_cp(T4, nasa_O2_a1[0], nasa_O2_a1[1], nasa_O2_a1[2], nasa_O2_a1[3], nasa_O2_a1[4]);

    const Scalar nasa_cp_O2_2 = R_O2 * nasa_cp(T2, nasa_O2_a2[0], nasa_O2_a2[1], nasa_O2_a2[2], nasa_O2_a2[3], nasa_O2_a2[4]);

    return_flag_temp = test_cp( species_name, index, nasa_cp_O2_1, T4, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_cp( species_name, index, nasa_cp_O2_2, T2, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;
  }

  // Test OH cp
  {
    unsigned int index = 2;
    const Scalar cp_1 = R_OH*cea_cp( T1, cea_OH_a1[0], cea_OH_a1[1], cea_OH_a1[2], cea_OH_a1[3], cea_OH_a1[4], cea_OH_a1[5], cea_OH_a1[6]);

    const Scalar cp_2 = R_OH*cea_cp( T2, cea_OH_a2[0], cea_OH_a2[1], cea_OH_a2[2], cea_OH_a2[3], cea_OH_a2[4], cea_OH_a2[5], cea_OH_a2[6]);

    const Scalar cp_3 = R_OH*cea_cp( T3, cea_OH_a3[0], cea_OH_a3[1], cea_OH_a3[2], cea_OH_a3[3], cea_OH_a3[4], cea_OH_a3[5], cea_OH_a3[6]);
    
    const Antioch::Species species = chem_mixture.species_list()[index];
    const std::string species_name = chem_mixture.species_inverse_name_map().find(species)->second;

    int return_flag_temp = 0;
    return_flag_temp = test_cp( species_name, index, cp_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_cp( species_name, index, cp_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_cp( species_name, index, cp_3, T3, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

// NASA [300 - 1000] [1000 - 5000]
    const Scalar nasa_cp_OH_1 = R_OH * nasa_cp(T4, nasa_OH_a1[0], nasa_OH_a1[1], nasa_OH_a1[2], nasa_OH_a1[3], nasa_OH_a1[4]);

    const Scalar nasa_cp_OH_2 = R_OH * nasa_cp(T2, nasa_OH_a2[0], nasa_OH_a2[1], nasa_OH_a2[2], nasa_OH_a2[3], nasa_OH_a2[4]);

    return_flag_temp = test_cp( species_name, index, nasa_cp_OH_1, T4, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_cp( species_name, index, nasa_cp_OH_2, T2, nasa_thermo );
  }


  // Test O cp
  {
    unsigned int index = 3;
    const Scalar cp_1 = R_O*cea_cp( T1, cea_O_a1[0], cea_O_a1[1], cea_O_a1[2], cea_O_a1[3], cea_O_a1[4], cea_O_a1[5], cea_O_a1[6]);

    const Scalar cp_2 = R_O*cea_cp( T2, cea_O_a2[0], cea_O_a2[1], cea_O_a2[2], cea_O_a2[3], cea_O_a2[4], cea_O_a2[5], cea_O_a2[6]);

    const Scalar cp_3 = R_O*cea_cp( T3, cea_O_a3[0], cea_O_a3[1], cea_O_a3[2], cea_O_a3[3], cea_O_a3[4], cea_O_a3[5], cea_O_a3[6]);
    
    const Antioch::Species species = chem_mixture.species_list()[index];
    const std::string species_name = chem_mixture.species_inverse_name_map().find(species)->second;

    int return_flag_temp = 0;
    return_flag_temp = test_cp( species_name, index, cp_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_cp( species_name, index, cp_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_cp( species_name, index, cp_3, T3, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

// NASA [300 - 1000] [1000 - 5000]
    const Scalar nasa_cp_O_1 = R_O * nasa_cp(T4, nasa_O_a1[0], nasa_O_a1[1], nasa_O_a1[2], nasa_O_a1[3], nasa_O_a1[4]);

    const Scalar nasa_cp_O_2 = R_O * nasa_cp(T2, nasa_O_a2[0], nasa_O_a2[1], nasa_O_a2[2], nasa_O_a2[3], nasa_O_a2[4]);

    return_flag_temp = test_cp( species_name, index, nasa_cp_O_1, T4, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_cp( species_name, index, nasa_cp_O_2, T2, nasa_thermo );
  }


  // Test H2 cp
  {
    unsigned int index = 4;
    const Scalar cp_1 = R_H2*cea_cp( T1, cea_H2_a1[0], cea_H2_a1[1], cea_H2_a1[2], cea_H2_a1[3], cea_H2_a1[4], cea_H2_a1[5], cea_H2_a1[6]);

    const Scalar cp_2 = R_H2*cea_cp( T2, cea_H2_a2[0], cea_H2_a2[1], cea_H2_a2[2], cea_H2_a2[3], cea_H2_a2[4], cea_H2_a2[5], cea_H2_a2[6]);

    const Scalar cp_3 = R_H2*cea_cp( T3, cea_H2_a3[0], cea_H2_a3[1], cea_H2_a3[2], cea_H2_a3[3], cea_H2_a3[4], cea_H2_a3[5], cea_H2_a3[6]);
    
    const Antioch::Species species = chem_mixture.species_list()[index];
    const std::string species_name = chem_mixture.species_inverse_name_map().find(species)->second;

    int return_flag_temp = 0;
    return_flag_temp = test_cp( species_name, index, cp_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_cp( species_name, index, cp_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_cp( species_name, index, cp_3, T3, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;


// NASA [300 - 1000] [1000 - 5000]
    const Scalar nasa_cp_H2_1 = R_H2 * nasa_cp(T4, nasa_H2_a1[0], nasa_H2_a1[1], nasa_H2_a1[2], nasa_H2_a1[3], nasa_H2_a1[4]);

    const Scalar nasa_cp_H2_2 = R_H2 * nasa_cp(T2, nasa_H2_a2[0], nasa_H2_a2[1], nasa_H2_a2[2], nasa_H2_a2[3], nasa_H2_a2[4]);

    return_flag_temp = test_cp( species_name, index, nasa_cp_H2_1, T4, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_cp( species_name, index, nasa_cp_H2_2, T2, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;
  }

// h, s, && g


  // Test N2 h
  {
    unsigned int index = 0;
    const Scalar h_N2_1 = R_N2 * T1 * cea_h( T1, cea_N2_a1[0], cea_N2_a1[1], cea_N2_a1[2], cea_N2_a1[3], cea_N2_a1[4], cea_N2_a1[5], cea_N2_a1[6],cea_N2_a1[8]);
    const Scalar s_N2_1 =             cea_s( T1, cea_N2_a1[0], cea_N2_a1[1], cea_N2_a1[2], cea_N2_a1[3], cea_N2_a1[4], cea_N2_a1[5], cea_N2_a1[6],cea_N2_a1[9]);
    const Scalar g_N2_1 = h_N2_1/(R_N2 * T1) - s_N2_1; // h/(R*T) - s/R

    const Scalar h_N2_2 = R_N2 * T2 * cea_h( T2, cea_N2_a2[0], cea_N2_a2[1], cea_N2_a2[2], cea_N2_a2[3], cea_N2_a2[4], cea_N2_a2[5], cea_N2_a2[6],cea_N2_a2[8]);
    const Scalar s_N2_2 =             cea_s( T2, cea_N2_a2[0], cea_N2_a2[1], cea_N2_a2[2], cea_N2_a2[3], cea_N2_a2[4], cea_N2_a2[5], cea_N2_a2[6],cea_N2_a2[9]);
    const Scalar g_N2_2 = h_N2_2/(R_N2 * T2) - s_N2_2;

    const Scalar h_N2_3 = R_N2 * T3 * cea_h( T3, cea_N2_a3[0], cea_N2_a3[1], cea_N2_a3[2], cea_N2_a3[3], cea_N2_a3[4], cea_N2_a3[5], cea_N2_a3[6],cea_N2_a3[8]);
    const Scalar s_N2_3 =             cea_s( T3, cea_N2_a3[0], cea_N2_a3[1], cea_N2_a3[2], cea_N2_a3[3], cea_N2_a3[4], cea_N2_a3[5], cea_N2_a3[6],cea_N2_a3[9]);
    const Scalar g_N2_3 = h_N2_3 /(R_N2 * T3) - s_N2_3;

    
    const Antioch::Species species = chem_mixture.species_list()[index];
    const std::string species_name = chem_mixture.species_inverse_name_map().find(species)->second;

    int return_flag_temp = 0;
    return_flag_temp = test_h( species_name, index, h_N2_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_h( species_name, index, h_N2_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_h( species_name, index, h_N2_3, T3, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_s( species_name, index, s_N2_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_s( species_name, index, s_N2_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_s( species_name, index, s_N2_3, T3, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_g( species_name, index, g_N2_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_g( species_name, index, g_N2_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_g( species_name, index, g_N2_3, T3, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

// NASA [300 - 1000] [1000 - 5000]
    const Scalar nasa_h_N2_1 = R_N2 * T4 * nasa_h(T4, nasa_N2_a1[0], nasa_N2_a1[1], nasa_N2_a1[2], nasa_N2_a1[3], nasa_N2_a1[4], nasa_N2_a1[5]);
    const Scalar nasa_s_N2_1 =             nasa_s(T4, nasa_N2_a1[0], nasa_N2_a1[1], nasa_N2_a1[2], nasa_N2_a1[3], nasa_N2_a1[4], nasa_N2_a1[6]);
    const Scalar nasa_g_N2_1 = nasa_h_N2_1/(R_N2 * T4) - nasa_s_N2_1;

    const Scalar nasa_h_N2_2 = R_N2 * T2 * nasa_h(T2, nasa_N2_a2[0], nasa_N2_a2[1], nasa_N2_a2[2], nasa_N2_a2[3], nasa_N2_a2[4], nasa_N2_a2[5]);
    const Scalar nasa_s_N2_2 =             nasa_s(T2, nasa_N2_a2[0], nasa_N2_a2[1], nasa_N2_a2[2], nasa_N2_a2[3], nasa_N2_a2[4], nasa_N2_a2[6]);
    const Scalar nasa_g_N2_2 = nasa_h_N2_2/(R_N2 * T2) - nasa_s_N2_2;

    return_flag_temp = test_h( species_name, index, nasa_h_N2_1, T4, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_h( species_name, index, nasa_h_N2_2, T2, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_s( species_name, index, nasa_s_N2_1, T4, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_s( species_name, index, nasa_s_N2_2, T2, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_g( species_name, index, nasa_g_N2_1, T4, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_g( species_name, index, nasa_g_N2_2, T2, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

  }

  // Test O2 h
  {
    unsigned int index = 1;
    const Scalar h_1 = R_O2*T1*cea_h( T1, cea_O2_a1[0], cea_O2_a1[1], cea_O2_a1[2], cea_O2_a1[3], cea_O2_a1[4], cea_O2_a1[5], cea_O2_a1[6],cea_O2_a1[8]);
    const Scalar s_1 =         cea_s( T1, cea_O2_a1[0], cea_O2_a1[1], cea_O2_a1[2], cea_O2_a1[3], cea_O2_a1[4], cea_O2_a1[5], cea_O2_a1[6],cea_O2_a1[9]);
    const Scalar g_1 = h_1/(R_O2*T1) - s_1;

    const Scalar h_2 = R_O2*T2*cea_h( T2, cea_O2_a2[0], cea_O2_a2[1], cea_O2_a2[2], cea_O2_a2[3], cea_O2_a2[4], cea_O2_a2[5], cea_O2_a2[6],cea_O2_a2[8]);
    const Scalar s_2 =         cea_s( T2, cea_O2_a2[0], cea_O2_a2[1], cea_O2_a2[2], cea_O2_a2[3], cea_O2_a2[4], cea_O2_a2[5], cea_O2_a2[6],cea_O2_a2[9]);
    const Scalar g_2 = h_2/(R_O2*T2) - s_2;

    const Scalar h_3 = R_O2*T3*cea_h( T3, cea_O2_a3[0], cea_O2_a3[1], cea_O2_a3[2], cea_O2_a3[3], cea_O2_a3[4], cea_O2_a3[5], cea_O2_a3[6],cea_O2_a3[8]);
    const Scalar s_3 =         cea_s( T3, cea_O2_a3[0], cea_O2_a3[1], cea_O2_a3[2], cea_O2_a3[3], cea_O2_a3[4], cea_O2_a3[5], cea_O2_a3[6],cea_O2_a3[9]);
    const Scalar g_3 = h_3/(R_O2*T3) - s_3;
    
    const Antioch::Species species = chem_mixture.species_list()[index];
    const std::string species_name = chem_mixture.species_inverse_name_map().find(species)->second;

    int return_flag_temp = 0;
    return_flag_temp = test_h( species_name, index, h_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_h( species_name, index, h_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_h( species_name, index, h_3, T3, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_s( species_name, index, s_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_s( species_name, index, s_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_s( species_name, index, s_3, T3, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_g( species_name, index, g_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_g( species_name, index, g_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_g( species_name, index, g_3, T3, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

// NASA [300 - 1000] [1000 - 5000]
    const Scalar nasa_h_O2_1 = R_O2 * T4 * nasa_h(T4, nasa_O2_a1[0], nasa_O2_a1[1], nasa_O2_a1[2], nasa_O2_a1[3], nasa_O2_a1[4], nasa_O2_a1[5]);
    const Scalar nasa_s_O2_1 =             nasa_s(T4, nasa_O2_a1[0], nasa_O2_a1[1], nasa_O2_a1[2], nasa_O2_a1[3], nasa_O2_a1[4], nasa_O2_a1[6]);
    const Scalar nasa_g_O2_1 = nasa_h_O2_1/(R_O2*T4) - nasa_s_O2_1;

    const Scalar nasa_h_O2_2 = R_O2 * T2 * nasa_h(T2, nasa_O2_a2[0], nasa_O2_a2[1], nasa_O2_a2[2], nasa_O2_a2[3], nasa_O2_a2[4], nasa_O2_a2[5]);
    const Scalar nasa_s_O2_2 =             nasa_s(T2, nasa_O2_a2[0], nasa_O2_a2[1], nasa_O2_a2[2], nasa_O2_a2[3], nasa_O2_a2[4], nasa_O2_a2[6]);
    const Scalar nasa_g_O2_2 = nasa_h_O2_2/(R_O2*T2) - nasa_s_O2_2;

    return_flag_temp = test_h( species_name, index, nasa_h_O2_1, T4, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_h( species_name, index, nasa_h_O2_2, T2, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_s( species_name, index, nasa_s_O2_1, T4, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_s( species_name, index, nasa_s_O2_2, T2, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_g( species_name, index, nasa_g_O2_1, T4, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_g( species_name, index, nasa_g_O2_2, T2, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

  }

  // Test OH h
  {
    unsigned int index = 2;
    const Scalar h_1 = R_OH*T1*cea_h( T1, cea_OH_a1[0], cea_OH_a1[1], cea_OH_a1[2], cea_OH_a1[3], cea_OH_a1[4], cea_OH_a1[5], cea_OH_a1[6],cea_OH_a1[8]);
    const Scalar s_1 =         cea_s( T1, cea_OH_a1[0], cea_OH_a1[1], cea_OH_a1[2], cea_OH_a1[3], cea_OH_a1[4], cea_OH_a1[5], cea_OH_a1[6],cea_OH_a1[9]);
    const Scalar g_1 = h_1/(R_OH*T1) - s_1;

    const Scalar h_2 = R_OH*T2*cea_h( T2, cea_OH_a2[0], cea_OH_a2[1], cea_OH_a2[2], cea_OH_a2[3], cea_OH_a2[4], cea_OH_a2[5], cea_OH_a2[6],cea_OH_a2[8]);
    const Scalar s_2 =         cea_s( T2, cea_OH_a2[0], cea_OH_a2[1], cea_OH_a2[2], cea_OH_a2[3], cea_OH_a2[4], cea_OH_a2[5], cea_OH_a2[6],cea_OH_a2[9]);
    const Scalar g_2 = h_2/(R_OH*T2) - s_2;

    const Scalar h_3 = R_OH*T3*cea_h( T3, cea_OH_a3[0], cea_OH_a3[1], cea_OH_a3[2], cea_OH_a3[3], cea_OH_a3[4], cea_OH_a3[5], cea_OH_a3[6],cea_OH_a3[8]);
    const Scalar s_3 =         cea_s( T3, cea_OH_a3[0], cea_OH_a3[1], cea_OH_a3[2], cea_OH_a3[3], cea_OH_a3[4], cea_OH_a3[5], cea_OH_a3[6],cea_OH_a3[9]);
    const Scalar g_3 = h_3/(R_OH*T3) - s_3;
    
    const Antioch::Species species = chem_mixture.species_list()[index];
    const std::string species_name = chem_mixture.species_inverse_name_map().find(species)->second;

    int return_flag_temp = 0;
    return_flag_temp = test_h( species_name, index, h_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_h( species_name, index, h_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_h( species_name, index, h_3, T3, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_s( species_name, index, s_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_s( species_name, index, s_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_s( species_name, index, s_3, T3, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_g( species_name, index, g_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_g( species_name, index, g_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_g( species_name, index, g_3, T3, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

// NASA [300 - 1000] [1000 - 5000]
    const Scalar nasa_h_OH_1 = R_OH * T4 * nasa_h(T4, nasa_OH_a1[0], nasa_OH_a1[1], nasa_OH_a1[2], nasa_OH_a1[3], nasa_OH_a1[4], nasa_OH_a1[5]);
    const Scalar nasa_s_OH_1 =             nasa_s(T4, nasa_OH_a1[0], nasa_OH_a1[1], nasa_OH_a1[2], nasa_OH_a1[3], nasa_OH_a1[4], nasa_OH_a1[6]);
    const Scalar nasa_g_OH_1 = nasa_h_OH_1/(R_OH*T4) - nasa_s_OH_1;

    const Scalar nasa_h_OH_2 = R_OH * T2 * nasa_h(T2, nasa_OH_a2[0], nasa_OH_a2[1], nasa_OH_a2[2], nasa_OH_a2[3], nasa_OH_a2[4], nasa_OH_a2[5]);
    const Scalar nasa_s_OH_2 =             nasa_s(T2, nasa_OH_a2[0], nasa_OH_a2[1], nasa_OH_a2[2], nasa_OH_a2[3], nasa_OH_a2[4], nasa_OH_a2[6]);
    const Scalar nasa_g_OH_2 = nasa_h_OH_2/(R_OH*T2) - nasa_s_OH_2;

    return_flag_temp = test_h( species_name, index, nasa_h_OH_1, T4, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_h( species_name, index, nasa_h_OH_2, T2, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_s( species_name, index, nasa_s_OH_1, T4, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_s( species_name, index, nasa_s_OH_2, T2, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_g( species_name, index, nasa_g_OH_1, T4, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_g( species_name, index, nasa_g_OH_2, T2, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;
  }

  // Test O h
  {
    unsigned int index = 3;
    const Scalar h_1 = R_O*T1*cea_h( T1, cea_O_a1[0], cea_O_a1[1], cea_O_a1[2], cea_O_a1[3], cea_O_a1[4], cea_O_a1[5], cea_O_a1[6],cea_O_a1[8]);
    const Scalar s_1 =        cea_s( T1, cea_O_a1[0], cea_O_a1[1], cea_O_a1[2], cea_O_a1[3], cea_O_a1[4], cea_O_a1[5], cea_O_a1[6],cea_O_a1[9]);
    const Scalar g_1 = h_1/(R_O*T1) - s_1;

    const Scalar h_2 = R_O*T2*cea_h( T2, cea_O_a2[0], cea_O_a2[1], cea_O_a2[2], cea_O_a2[3], cea_O_a2[4], cea_O_a2[5], cea_O_a2[6],cea_O_a2[8]);
    const Scalar s_2 =        cea_s( T2, cea_O_a2[0], cea_O_a2[1], cea_O_a2[2], cea_O_a2[3], cea_O_a2[4], cea_O_a2[5], cea_O_a2[6],cea_O_a2[9]);
    const Scalar g_2 = h_2/(R_O*T2) - s_2;

    const Scalar h_3 = R_O*T3*cea_h( T3, cea_O_a3[0], cea_O_a3[1], cea_O_a3[2], cea_O_a3[3], cea_O_a3[4], cea_O_a3[5], cea_O_a3[6],cea_O_a3[8]);
    const Scalar s_3 =        cea_s( T3, cea_O_a3[0], cea_O_a3[1], cea_O_a3[2], cea_O_a3[3], cea_O_a3[4], cea_O_a3[5], cea_O_a3[6],cea_O_a3[9]);
    const Scalar g_3 = h_3/(R_O*T3) - s_3;
    
    const Antioch::Species species = chem_mixture.species_list()[index];
    const std::string species_name = chem_mixture.species_inverse_name_map().find(species)->second;

    int return_flag_temp = 0;
    return_flag_temp = test_h( species_name, index, h_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_h( species_name, index, h_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_h( species_name, index, h_3, T3, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_s( species_name, index, s_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_s( species_name, index, s_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_s( species_name, index, s_3, T3, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_g( species_name, index, g_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_g( species_name, index, g_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_g( species_name, index, g_3, T3, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

// NASA [300 - 1000] [1000 - 5000]
    const Scalar nasa_h_O_1 = R_O * T4 * nasa_h(T4, nasa_O_a1[0], nasa_O_a1[1], nasa_O_a1[2], nasa_O_a1[3], nasa_O_a1[4], nasa_O_a1[5]);
    const Scalar nasa_s_O_1 =            nasa_s(T4, nasa_O_a1[0], nasa_O_a1[1], nasa_O_a1[2], nasa_O_a1[3], nasa_O_a1[4], nasa_O_a1[6]);
    const Scalar nasa_g_O_1 = nasa_h_O_1/(R_O*T4) - nasa_s_O_1;

    const Scalar nasa_h_O_2 = R_O * T2 * nasa_h(T2, nasa_O_a2[0], nasa_O_a2[1], nasa_O_a2[2], nasa_O_a2[3], nasa_O_a2[4], nasa_O_a2[5]);
    const Scalar nasa_s_O_2 =            nasa_s(T2, nasa_O_a2[0], nasa_O_a2[1], nasa_O_a2[2], nasa_O_a2[3], nasa_O_a2[4], nasa_O_a2[6]);
    const Scalar nasa_g_O_2 = nasa_h_O_2/(R_O*T2) - nasa_s_O_2;

    return_flag_temp = test_h( species_name, index, nasa_h_O_1, T4, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_h( species_name, index, nasa_h_O_2, T2, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_s( species_name, index, nasa_s_O_1, T4, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_s( species_name, index, nasa_s_O_2, T2, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_g( species_name, index, nasa_g_O_1, T4, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_g( species_name, index, nasa_g_O_2, T2, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;
  }

  // Test H2 h
  {
    unsigned int index = 4;
    const Scalar h_1 = R_H2*T1*cea_h( T1, cea_H2_a1[0], cea_H2_a1[1], cea_H2_a1[2], cea_H2_a1[3], cea_H2_a1[4], cea_H2_a1[5], cea_H2_a1[6],cea_H2_a1[8]);
    const Scalar s_1 =         cea_s( T1, cea_H2_a1[0], cea_H2_a1[1], cea_H2_a1[2], cea_H2_a1[3], cea_H2_a1[4], cea_H2_a1[5], cea_H2_a1[6],cea_H2_a1[9]);
    const Scalar g_1 = h_1/(R_H2*T1) - s_1;

    const Scalar h_2 = R_H2*T2*cea_h( T2, cea_H2_a2[0], cea_H2_a2[1], cea_H2_a2[2], cea_H2_a2[3], cea_H2_a2[4], cea_H2_a2[5], cea_H2_a2[6],cea_H2_a2[8]);
    const Scalar s_2 =         cea_s( T2, cea_H2_a2[0], cea_H2_a2[1], cea_H2_a2[2], cea_H2_a2[3], cea_H2_a2[4], cea_H2_a2[5], cea_H2_a2[6],cea_H2_a2[9]);
    const Scalar g_2 = h_2/(R_H2*T2) - s_2;

    const Scalar h_3 = R_H2*T3*cea_h( T3, cea_H2_a3[0], cea_H2_a3[1], cea_H2_a3[2], cea_H2_a3[3], cea_H2_a3[4], cea_H2_a3[5], cea_H2_a3[6],cea_H2_a3[8]);
    const Scalar s_3 =         cea_s( T3, cea_H2_a3[0], cea_H2_a3[1], cea_H2_a3[2], cea_H2_a3[3], cea_H2_a3[4], cea_H2_a3[5], cea_H2_a3[6],cea_H2_a3[9]);
    const Scalar g_3 = h_3/(R_H2*T3) - s_3;
    
    const Antioch::Species species = chem_mixture.species_list()[index];
    const std::string species_name = chem_mixture.species_inverse_name_map().find(species)->second;

    int return_flag_temp = 0;
    return_flag_temp = test_h( species_name, index, h_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_h( species_name, index, h_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_h( species_name, index, h_3, T3, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_s( species_name, index, s_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_s( species_name, index, s_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_s( species_name, index, s_3, T3, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_g( species_name, index, g_1, T1, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_g( species_name, index, g_2, T2, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_g( species_name, index, g_3, T3, thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

return return_flag;

// NASA [300 - 1000] [1000 - 5000]
    const Scalar nasa_h_H2_1 = R_H2 * T4 * nasa_h(T4, nasa_H2_a1[0], nasa_H2_a1[1], nasa_H2_a1[2], nasa_H2_a1[3], nasa_H2_a1[4], nasa_H2_a1[5]);
    const Scalar nasa_s_H2_1 =             nasa_s(T4, nasa_H2_a1[0], nasa_H2_a1[1], nasa_H2_a1[2], nasa_H2_a1[3], nasa_H2_a1[4], nasa_H2_a1[6]);
    const Scalar nasa_g_H2_1 = nasa_h_H2_1/(R_H2 * T4) - nasa_s_H2_1;

    const Scalar nasa_h_H2_2 = R_H2 * T2 * nasa_h(T2, nasa_H2_a2[0], nasa_H2_a2[1], nasa_H2_a2[2], nasa_H2_a2[3], nasa_H2_a2[4], nasa_H2_a2[5]);
    const Scalar nasa_s_H2_2 =             nasa_s(T2, nasa_H2_a2[0], nasa_H2_a2[1], nasa_H2_a2[2], nasa_H2_a2[3], nasa_H2_a2[4], nasa_H2_a2[6]);
    const Scalar nasa_g_H2_2 = nasa_h_H2_2/(R_H2 * T2) - nasa_s_H2_2;

    return_flag_temp = test_h( species_name, index, nasa_h_H2_1, T4, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_h( species_name, index, nasa_h_H2_2, T2, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_s( species_name, index, nasa_s_H2_1, T4, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_s( species_name, index, nasa_s_H2_2, T2, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_g( species_name, index, nasa_g_H2_1, T4, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;

    return_flag_temp = test_g( species_name, index, nasa_g_H2_2, T2, nasa_thermo );
    if( return_flag_temp != 0 ) return_flag = 1;
  }

  return return_flag;
}


int main(int argc, char* argv[])
{
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify NASA thermodynamics set input file." << std::endl;
      antioch_error();
    }
// We're not getting the full long double precision yet?
  return (tester<double>(std::string(argv[1])) ||
//          tester<long double>(std::string(argv[1])) ||
          tester<float>(std::string(argv[1])));
}
