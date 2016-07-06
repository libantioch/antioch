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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// C++
#include <limits>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>

// Antioch
#include "antioch/vector_utils.h"

#include "antioch/antioch_asserts.h"
#include "antioch/chemical_species.h"
#include "antioch/chemical_mixture.h"
#include "antioch/reaction_set.h"
#include "antioch/read_reaction_set_data.h"
#include "antioch/nasa_mixture_parsing.h"
#include "antioch/nasa_evaluator.h"
#include "antioch/kinetics_evaluator.h"


template <typename Scalar>
std::vector<Scalar> temp_thermo()
{
  std::vector<Scalar> Ts(3,0);
  Ts[0] = 200;
  Ts[1] = 1000;
  Ts[2] = 6000;
  return Ts;
}

template <typename Scalar>
std::vector<Scalar> H_thermo(const Scalar & T)
{
  std::vector<std::vector<Scalar> > out(3,std::vector<Scalar>(10,0));
  out[0][0] =  0.00000000e+00;  
  out[0][1] =  0.00000000e+00;
  out[0][2] =  2.50000000e+00;
  out[0][3] =  0.00000000e+00;
  out[0][4] =  0.00000000e+00;
  out[0][5] =  0.00000000e+00;
  out[0][6] =  0.00000000e+00;
  out[0][7] =  0.00000000e+00;
  out[0][8] =  2.54737080e+04;
  out[0][9] = -4.46682853e-01;

  out[1][0] =  6.07877425e+01;
  out[1][1] = -1.81935442e-01;
  out[1][2] =  2.50021182e+00;
  out[1][3] = -1.22651286e-07;
  out[1][4] =  3.73287633e-11;
  out[1][5] = -5.68774456e-15;
  out[1][6] =  3.41021020e-19;
  out[1][7] =  0.00000000e+00;
  out[1][8] =  2.54748640e+04;
  out[1][9] = -4.48191777e-01;

  out[2][0] =  2.17375769e+08;
  out[2][1] = -1.31203540e+05;
  out[2][2] =  3.39917420e+01;
  out[2][3] = -3.81399968e-03;
  out[2][4] =  2.43285484e-07;
  out[2][5] = -7.69427554e-12;
  out[2][6] =  9.64410563e-17;
  out[2][7] =  0.00000000e+00;
  out[2][8] =  1.06763809e+06;
  out[2][9] = -2.74230105e+02;

  return (T <= temp_thermo<Scalar>()[1])?out[0]:(T <= temp_thermo<Scalar>()[2])?out[1]:out[2];
}

template <typename Scalar>
std::vector<Scalar> H2_thermo(const Scalar & T)
{
  std::vector<std::vector<Scalar> > out(3,std::vector<Scalar>(10,0));

  out[0][0] =  4.07832281e+04;
  out[0][1] = -8.00918545e+02;
  out[0][2] =  8.21470167e+00;
  out[0][3] = -1.26971436e-02;
  out[0][4] =  1.75360493e-05;
  out[0][5] = -1.20286016e-08;
  out[0][6] =  3.36809316e-12;
  out[0][7] =  0.00000000e+00;
  out[0][8] =  2.68248438e+03;
  out[0][9] = -3.04378866e+01;

  out[1][0] =  5.60812338e+05;
  out[1][1] = -8.37149134e+02;
  out[1][2] =  2.97536304e+00;
  out[1][3] =  1.25224993e-03;
  out[1][4] = -3.74071842e-07;
  out[1][5] =  5.93662825e-11;
  out[1][6] = -3.60699573e-15;
  out[1][7] =  0.00000000e+00;
  out[1][8] =  5.33981585e+03;
  out[1][9] = -2.20276405e+00;

  out[2][0] = 4.96671613e+08;
  out[2][1] = -3.14744812e+05;
  out[2][2] = 7.98388750e+01;
  out[2][3] = -8.41450419e-03;
  out[2][4] = 4.75306044e-07;
  out[2][5] = -1.37180973e-11;
  out[2][6] = 1.60537460e-16;
  out[2][7] = 0.00000000e+00;
  out[2][8] = 2.48835466e+06;
  out[2][9] = -6.69552419e+02;

  return (T <= temp_thermo<Scalar>()[1])?out[0]:(T <= temp_thermo<Scalar>()[2])?out[1]:out[2];
}

template <typename Scalar>
std::vector<Scalar> H2O2_thermo(const Scalar & T)
{
  std::vector<std::vector<Scalar> > out(2,std::vector<Scalar>(10,0));

  out[0][0] = -9.279533580E+04;
  out[0][1] =  1.564748385E+03;
  out[0][2] = -5.976460140E+00;
  out[0][3] =  3.270744520E-02;
  out[0][4] = -3.932193260E-05;
  out[0][5] =  2.509255235E-08;
  out[0][6] = -6.465045290E-12;
  out[0][7] =  0.000000000E+00;
  out[0][8] = -2.494004728E+04;
  out[0][9] =  5.877174180E+01;

  out[1][0] =  1.489428027E+06;
  out[1][1] = -5.170821780E+03;
  out[1][2] =  1.128204970E+01;
  out[1][3] = -8.042397790E-05;
  out[1][4] = -1.818383769E-08;
  out[1][5] =  6.947265590E-12;
  out[1][6] = -4.827831900E-16;
  out[1][7] =  0.000000000E+00;
  out[1][8] =  1.418251038E+04;
  out[1][9] = -4.650855660E+01;

  return (T <= temp_thermo<Scalar>()[1])?out[0]:out[1];
}
										 
template <typename Scalar>
std::vector<Scalar> HO2_thermo(const Scalar & T)
{
  std::vector<std::vector<Scalar> > out(2,std::vector<Scalar>(10,0));

  out[0][0] = -7.598882540E+04;
  out[0][1] =  1.329383918E+03;
  out[0][2] = -4.677388240E+00;
  out[0][3] =  2.508308202E-02;
  out[0][4] = -3.006551588E-05;
  out[0][5] =  1.895600056E-08;
  out[0][6] = -4.828567390E-12;
  out[0][7] =  0.000000000E+00;
  out[0][8] = -5.873350960E+03;
  out[0][9] =  5.193602140E+01;

  out[1][0] = -1.810669724E+06;
  out[1][1] =  4.963192030E+03;
  out[1][2] = -1.039498992E+00;
  out[1][3] =  4.560148530E-03;
  out[1][4] = -1.061859447E-06;
  out[1][5] =  1.144567878E-10;
  out[1][6] = -4.763064160E-15;
  out[1][7] =  0.000000000E+00;
  out[1][8] = -3.200817190E+04;
  out[1][9] =  4.066850920E+01;

  return (T <= temp_thermo<Scalar>()[1])?out[0]:out[1];
}

template <typename Scalar>
Scalar g(const std::vector<Scalar> & thermo, const Scalar & T)
{
  return - thermo[0] / ( 2 * T * T)  
         + thermo[1] * ( 1 + Antioch::ant_log(T) ) / T  
         + thermo[2] * ( 1 - Antioch::ant_log(T) )
         - thermo[3] * T / 2
         - thermo[4] * T * T / 6
         - thermo[5] * T * T * T / 12 
         - thermo[6] * T * T * T * T / 20
         + thermo[8] / T
         - thermo[9];
}

template <typename Scalar>
Scalar G(const Scalar & T)
{
  return g(H2O2_thermo(T),T) + g(H_thermo(T),T) - g(HO2_thermo(T),T) - g(H2_thermo(T),T);
}

template <typename Scalar>
int check_test(const Scalar &exact,const Scalar &cal,const std::string &words)
{
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 250;

  if(std::abs(exact - cal)/exact > tol)
  {
     std::cout << std::scientific << std::setprecision(20)
               << "Erreur in tests of "  << words                       << "\n"
               << "Calculated value is " << cal                         << "\n"
               << "Exact value is "      << exact                       << "\n"
               << "Relative error is "   << std::abs(exact - cal)/exact << "\n"
               << "Tolerance is "        << tol << std::endl;
    return 1;
  }
  return 0;
}

template <typename Scalar>
int checker(const Scalar & net_rates_exact, 
            const Scalar & kfwd_const_exact,  const Scalar & kfwd_exact, const Scalar & fwd_conc_exact, 
            const Scalar & kbkwd_const_exact, const Scalar & kbkwd_exact, const Scalar & bkwd_conc_exact,
            const Scalar & net_rates, 
            const Scalar & kfwd_const,  const Scalar & kfwd,  const Scalar & fwd_conc, 
            const Scalar & kbkwd_const, const Scalar & kbkwd, const Scalar & bkwd_conc,
            const Scalar & Temp)
{

  int return_flag(0);

  std::stringstream os;
  os << Temp << "K";

  return_flag = check_test(fwd_conc_exact,fwd_conc,"concentrations forward at " + os.str())       ||
                return_flag;

  return_flag = check_test(kfwd_const_exact,kfwd_const,"rate constant forward at " + os.str())    ||
                return_flag;

  return_flag = check_test(kfwd_exact,kfwd,"rate forward at " + os.str())                         ||
                return_flag;

  return_flag = check_test(bkwd_conc_exact,bkwd_conc,"concentrations backward at " + os.str())    ||
                return_flag;

  return_flag = check_test(kbkwd_const_exact,kbkwd_const,"rate constant backward at " + os.str()) ||
                return_flag;

  return_flag = check_test(kbkwd_exact,kbkwd,"rate backward at " + os.str())                      ||
                return_flag;

  return_flag = check_test(net_rates_exact,net_rates,"net rate at " + os.str())                   ||
                return_flag;

  return return_flag;
}

template <typename Scalar>
int test_type(const std::string& input_name, const Antioch::ParsingType & inputType)
{
  using std::abs;

  std::vector<std::string> species_str_list;
  const unsigned int n_species = 4;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "H" );
  species_str_list.push_back( "H2" );
  species_str_list.push_back( "H2O2" );
  species_str_list.push_back( "HO2" );

  // kinetics
  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );
  Antioch::ReactionSet<Scalar> reaction_set( chem_mixture );
  Antioch::read_reaction_set_data<Scalar>( input_name, true, reaction_set, inputType );

  // thermo
  Antioch::NASAThermoMixture<Scalar,Antioch::CEACurveFit<Scalar> > thermo_mixture(chem_mixture); 
  Antioch::read_nasa_mixture_data( thermo_mixture ); // default
  Antioch::NASAEvaluator<Scalar,Antioch::CEACurveFit<Scalar> > thermo(thermo_mixture); 

  const Scalar T = 1500.0; // K
  const Scalar P = 1.0e5; // Pa

  // Mass fractions
  std::vector<Scalar> Y(n_species,0.25);

  const Scalar R_mix = chem_mixture.R(Y); // get R_tot in J.kg-1.K-1

  std::vector<Scalar> molar_densities(n_species,0);

  std::vector<Scalar> h_RT_minus_s_R(n_species);

  Scalar net_rates_exact(0),
         kfwd_const_exact(0),
         kbkwd_const_exact(0),
         kfwd_exact(0),
         kbkwd_exact(0),
         fwd_conc_exact(0),
         bkwd_conc_exact(0);

  std::vector<Scalar> net_rates,kfwd_const,kbkwd_const,kfwd,kbkwd,fwd_conc,bkwd_conc;

  // golden values from bc
  kfwd_const_exact  = 3341803.012517167957974472239343341385324934527797794;
  fwd_conc_exact    = 0.13473900180907885567483372493084259843941997896096;
  kfwd_exact        = 450271.20214913586326479266723410049955538508784467327734;
  bkwd_conc_exact   = 0.06329787652743180276417606259788367794531198101580;
  kbkwd_const_exact = 4327.11114439947889901291861278625576793065735934819023;
  kbkwd_exact       = 273.89694693867234146591036506821238311545088830214564;
  net_rates_exact   = 449997.30520219719092332675686903228717226963695635680885;


  Scalar rho = P/(R_mix*T); // kg.m-3
  chem_mixture.molar_densities(rho,Y,molar_densities);
  const Antioch::KineticsConditions<Scalar> conditions(T);
  const Antioch::TempCache<Scalar> Cache(T);
  thermo.h_RT_minus_s_R(Cache,h_RT_minus_s_R);
  reaction_set.get_reactive_scheme(conditions,molar_densities,h_RT_minus_s_R,net_rates,
                                   kfwd_const,kbkwd_const,kfwd,kbkwd,fwd_conc,bkwd_conc);

  int return_flag = checker(net_rates_exact, kfwd_const_exact, kfwd_exact, fwd_conc_exact, kbkwd_const_exact, kbkwd_exact, bkwd_conc_exact,
                            net_rates[0],    kfwd_const[0],    kfwd[0],    fwd_conc[0],    kbkwd_const[0],    kbkwd[0],    bkwd_conc[0], T);

  const Scalar Rcal = Antioch::Constants::R_universal<Scalar>() * Antioch::Constants::R_universal_unit<Scalar>().factor_to_some_unit("cal/mol/K");
  const Scalar fac(1e-6);

  for(Scalar Temp = 210; Temp < 5990; Temp += 10){

    rho = P/(R_mix*Temp); // kg.m-3

    molar_densities.resize(4,0);
    chem_mixture.molar_densities(rho,Y,molar_densities);

    const Antioch::KineticsConditions<Scalar> conditionsTemp(Temp);
    const Antioch::TempCache<Scalar> CacheTemp(Temp);
    h_RT_minus_s_R.resize(4,0);
    thermo.h_RT_minus_s_R(CacheTemp,h_RT_minus_s_R);

    net_rates.clear();
    kfwd_const.clear();
    kbkwd_const.clear();
    kfwd.clear();
    kbkwd.clear();
    fwd_conc.clear();
    bkwd_conc.clear();
    reaction_set.get_reactive_scheme(conditionsTemp,molar_densities,h_RT_minus_s_R,net_rates,
                                     kfwd_const,kbkwd_const,kfwd,kbkwd,fwd_conc,bkwd_conc);

// H2O2 + H <=> HO2 + H2
// k = Arrhenius(4.82e13 (cm3/mol/s), 7.95e3 (cal/mol) )
// m_H2O2 = 1.5
// m_H = 0.5
// m_HO2 = 2
// m_H2 = 1
    fwd_conc_exact    = Antioch::ant_pow(molar_densities[2],1.5) * Antioch::ant_pow(molar_densities[0],0.5);
    bkwd_conc_exact   = Antioch::ant_pow(molar_densities[3],2) * molar_densities[1];
    kfwd_const_exact  = 4.82e13 * fac * std::exp(-7.95e3/(Rcal * Temp)); // Arrhenius
    kfwd_exact        = kfwd_const_exact * fwd_conc_exact;
    kbkwd_const_exact = kfwd_const_exact / Antioch::ant_exp(G(Temp));
    kbkwd_exact       = kbkwd_const_exact * bkwd_conc_exact;
    net_rates_exact   = kfwd_exact - kbkwd_exact;

    return_flag = checker(net_rates_exact, kfwd_const_exact, kfwd_exact, fwd_conc_exact, kbkwd_const_exact, kbkwd_exact, bkwd_conc_exact,
                          net_rates[0],    kfwd_const[0],    kfwd[0],    fwd_conc[0],    kbkwd_const[0],    kbkwd[0],    bkwd_conc[0], Temp) ||
                  return_flag;
  }

  return return_flag;
}

template <typename Scalar>
int tester(const std::string& input_name_xml, const std::string& input_name_ck)
{
  return test_type<Scalar>(input_name_xml, Antioch::ParsingType::XML) ||
    test_type<Scalar>(input_name_ck, Antioch::ParsingType::CHEMKIN);
}

int main(int argc, char* argv[])
{
  // Check command line count.
  if( argc < 3 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify reaction set XML input file and reaction set ChemKin input file." << std::endl;
      antioch_error();
    }

  return (tester<float>(std::string(argv[1]), std::string(argv[2]) ) ||
          tester<double>(std::string(argv[1]), std::string(argv[2]) )  /*||
          tester<long double>(std::string(argv[1]), std::string(argv[2]) )*/
         ); 
}
