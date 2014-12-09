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
#include <limits>
#include <iomanip>
#include <string>
#include <vector>

// Antioch
#include "antioch/vector_utils.h"

#include "antioch/antioch_asserts.h"
#include "antioch/chemical_species.h"
#include "antioch/chemical_mixture.h"
#include "antioch/reaction_set.h"
#include "antioch/read_reaction_set_data.h"
#include "antioch/cea_thermo.h"
#include "antioch/kinetics_evaluator.h"


template <typename Scalar>
int check_test(const Scalar &exact,const Scalar &cal,const std::string &words, unsigned int rxn)
{
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;

  if(std::abs(exact - cal)/exact > tol)
  {
     std::cout << std::scientific << std::setprecision(20)
               << "Erreur in tests of "  << words                       << "\n"
               << "of reaction #"        << rxn                         << "\n"
               << "Calculated value is " << cal                         << "\n"
               << "Exact value is "      << exact                       << "\n"
               << "Relative error is "   << std::abs(exact - cal)/exact << "\n"
               << "Tolerance is "        << tol << std::endl;
    return 1;
  }
  return 0;
}

template <typename Scalar>
int tester(const std::string& input_name)
{
  using std::abs;

  std::vector<std::string> species_str_list;
  const unsigned int n_species = 5;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "O2" );
  species_str_list.push_back( "N" );
  species_str_list.push_back( "O" );
  species_str_list.push_back( "NO" );

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );
  Antioch::ReactionSet<Scalar> reaction_set( chem_mixture );
  Antioch::CEAThermodynamics<Scalar> thermo( chem_mixture );

  Antioch::read_reaction_set_data_xml<Scalar>( input_name, true, reaction_set );

  const Scalar T = 1500.0; // K
  const Scalar P = 1.0e5; // Pa

  const Antioch::KineticsConditions<Scalar> conditions(T);

  // Mass fractions
  std::vector<Scalar> Y(n_species,0.2);

  const Scalar R_mix = chem_mixture.R(Y); // get R_tot in J.kg-1.K-1

  const Scalar rho = P/(R_mix*T); // kg.m-3

  std::vector<Scalar> molar_densities(n_species,0.0);
  chem_mixture.molar_densities(rho,Y,molar_densities);

  std::vector<Scalar> h_RT_minus_s_R(n_species);

  typedef typename Antioch::CEAThermodynamics<Scalar>::template Cache<Scalar> Cache;
  thermo.h_RT_minus_s_R(Cache(T),h_RT_minus_s_R);

  std::vector<Scalar> net_rates,kfwd_const,kbkwd_const,kfwd,kbkwd,fwd_conc,bkwd_conc;

  reaction_set.get_reactive_scheme(conditions,molar_densities,h_RT_minus_s_R,net_rates,
                                   kfwd_const,kbkwd_const,kfwd,kbkwd,fwd_conc,bkwd_conc);

  std::vector<Scalar> net_rates_exact,
                      kfwd_const_exact,
                      kbkwd_const_exact,
                      kfwd_exact,
                      kbkwd_exact,
                      fwd_conc_exact,
                      bkwd_conc_exact;
  net_rates_exact.resize(5,0.L);
  kfwd_const_exact.resize(5,0.L);
  kbkwd_const_exact.resize(5,0.L);
  kfwd_exact.resize(5,0.L);
  kbkwd_exact.resize(5,0.L);
  fwd_conc_exact.resize(5,1.L);
  bkwd_conc_exact.resize(5,1.L);

  const Scalar Rcal = Antioch::Constants::R_universal<Scalar>() * Antioch::Constants::R_universal_unit<Scalar>().factor_to_some_unit("cal/mol/K");
  const Scalar P0_RT(1.0e5/Antioch::Constants::R_universal<Scalar>()/T);
//N2 + M <=> 2 N + M
  kfwd_const_exact[0] = 7e15L * std::pow(T,-1.6L) * std::exp(-224801.3L/(Rcal * T)) * //Kooij
                        (4.2857L * (molar_densities[2] + molar_densities[3]) // N and O
                         + molar_densities[0] + molar_densities[1] + molar_densities[4] ); //N2, O2, NO
  fwd_conc_exact[0] = molar_densities[0]; //N2
  bkwd_conc_exact[0] = std::pow(molar_densities[2],2); //2 N

//O2 + M <=> 2 O + M
  kfwd_const_exact[1] = 2e15L * std::pow(T,-1.5L) * std::exp(-117881.7L/(Rcal * T)) * //Kooij
                        (5.0L * (molar_densities[2] + molar_densities[3]) // N and O
                         + molar_densities[0] + molar_densities[1] + molar_densities[4] ); //N2, O2, NO
  fwd_conc_exact[1] = molar_densities[1]; //O2
  bkwd_conc_exact[1] = std::pow(molar_densities[3],2); //2 O

//NO + M <=> N + O + M
  kfwd_const_exact[2] = 5e9L * std::exp(-149943.0L/(Rcal * T)) * //Kooij
                        (22.0L * (molar_densities[2] + molar_densities[3] + molar_densities[4]) // N, O and NO
                         + molar_densities[0] + molar_densities[1]); //N2, O2
  fwd_conc_exact[2] = molar_densities[4]; //NO
  bkwd_conc_exact[2] = molar_densities[2] * molar_densities[3]; //N + O

//N2 + O => NO + N
  kfwd_const_exact[3] = 5.7e9L * std::pow(T,0.42) * std::exp(-85269.6L/(Rcal * T)); //Kooij
  fwd_conc_exact[3] = molar_densities[0] * molar_densities[3]; //N2 + O

//NO + O => O2 + N
  kfwd_const_exact[4] = 8.4e9L * std::exp(-38526.0L/(Rcal * T)); //Kooij
  fwd_conc_exact[4] = molar_densities[4] * molar_densities[3]; //NO + O

  for(unsigned int rxn = 0; rxn < 5; rxn++)
  {
    if(rxn < 3)
    {
      kbkwd_const_exact[rxn] = kfwd_const_exact[rxn]/reaction_set.reaction(rxn).equilibrium_constant(P0_RT,h_RT_minus_s_R);
      kbkwd_exact[rxn] = kbkwd_const_exact[rxn] * bkwd_conc_exact[rxn];
    }
    kfwd_exact[rxn] = kfwd_const_exact[rxn] * fwd_conc_exact[rxn];
    net_rates_exact[rxn] = kfwd_exact[rxn] - kbkwd_exact[rxn];
  }

  int return_flag(0);
  for(unsigned int rxn = 0; rxn < 5; rxn++)
  {
     return_flag = check_test(net_rates_exact[rxn],net_rates[rxn],"net rates",rxn)                  ||
                   return_flag;
     return_flag = check_test(kfwd_const_exact[rxn],kfwd_const[rxn],"rate constant forward",rxn)    ||
                   return_flag;
     return_flag = check_test(kfwd_exact[rxn],kfwd[rxn],"rate forward",rxn)                         ||
                   return_flag;
     return_flag = check_test(fwd_conc_exact[rxn],fwd_conc[rxn],"concentrations forward",rxn)       ||
                   return_flag;

     if(rxn < 3)
     {
        return_flag = check_test(kbkwd_const_exact[rxn],kbkwd_const[rxn],"rate constant backward",rxn) ||
                      return_flag;
        return_flag = check_test(kbkwd_exact[rxn],kbkwd[rxn],"rate backward",rxn)                      ||
                      return_flag;
        return_flag = check_test(bkwd_conc_exact[rxn],bkwd_conc[rxn],"concentrations backward",rxn)    ||
                      return_flag;
     }
  }
  return return_flag;
}

int main(int argc, char* argv[])
{
  // Check command line count.
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify reaction set XML input file." << std::endl;
      antioch_error();
    }

  return (tester<float>(std::string(argv[1])) ||
          tester<double>(std::string(argv[1])));/* ||
          tester<long double>(std::string(argv[1]))
          ); */
}
