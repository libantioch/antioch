//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014 Paul T. Bauman, Sylvain Plessis, Roy H. Stonger
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

//Antioch
#include "antioch/vector_utils.h"
#include "antioch/reaction.h"
#include "antioch/reaction_enum.h"
#include "antioch/reaction_parsing.h"
#include "antioch/kinetics_parsing.h"
#include "antioch/chemical_species.h"
#include "antioch/chemical_mixture.h"
#include "antioch/physical_constants.h"
#include "antioch/cea_thermo.h"
#include "antioch/cmath_shims.h"

//C++
#include <iomanip>
#include <string>
#include <limits>

template <typename Scalar>
int tester(const std::string& testname)
{
  using std::abs;
  using std::exp;
  using std::pow;

  std::vector<std::string> species_str_list;
  const unsigned int n_species = 4;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "O" );
  species_str_list.push_back( "NO" );
  species_str_list.push_back( "N" );

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );
  Antioch::CEAThermodynamics<Scalar> thermo( chem_mixture );

 //fill on hand reaction_set
  const std::string equation("N2 + O [=] NO + N");
  const Scalar A = 5.7e+9L;
  const Scalar b = 0.42L;
  const Scalar E = 85269.6L;// cal/mol
  const Antioch::ReactionType::ReactionType typeReaction(Antioch::ReactionType::ELEMENTARY);
  const Antioch::KineticsModel::KineticsModel kineticsModel(Antioch::KineticsModel::KOOIJ);
  const bool reversible(true);

  const Scalar T = 1500.0L; // K
  const Scalar P = 1.0e5L; // Pa

  const Scalar P0_RT(1.0e5/Antioch::Constants::R_universal<Scalar>()/T);
  std::vector<Scalar> h_RT_minus_s_R(n_species);

  typedef typename Antioch::CEAThermodynamics<Scalar>::template Cache<Scalar> Cache;
  thermo.h_RT_minus_s_R(Cache(T),h_RT_minus_s_R);
  // Molar densities, equimolar
  std::vector<Scalar> molar_densities(n_species , P/Antioch::Constants::R_universal<Scalar>()/T/Scalar(n_species));
// Theory
  Scalar k         = A * Antioch::ant_pow(T,b) * Antioch::ant_exp(-E/(1.9858775L*T));
  Scalar Keq       = Antioch::ant_exp(h_RT_minus_s_R[0] + h_RT_minus_s_R[1] - h_RT_minus_s_R[2] - h_RT_minus_s_R[3]);
  Scalar kback     = k/Keq;
  Scalar rate_fwd  = k * molar_densities[0] * molar_densities[1];
  Scalar rate_back = kback * molar_densities[2] * molar_densities[3];
  Scalar net_rate  = rate_fwd - rate_back;

//reversible reaction  
  Antioch::Reaction<Scalar>* my_rxn = Antioch::build_reaction<Scalar>(n_species, equation,reversible,typeReaction,kineticsModel);
  std::vector<Scalar> data;
  data.push_back(A);
  data.push_back(b);
  data.push_back(E);
  data.push_back(1.L);
  data.push_back(1.9858775L);
  Antioch::KineticsType<Scalar>* rate = Antioch::build_rate<Scalar>(data,kineticsModel);
  my_rxn->add_forward_rate(rate);
  my_rxn->add_reactant("N2",chem_mixture.species_list_map().at(Antioch::Species::N2),1);
  my_rxn->add_reactant("O" ,chem_mixture.species_list_map().at(Antioch::Species::O) ,1);
  my_rxn->add_product ("NO",chem_mixture.species_list_map().at(Antioch::Species::NO),1);
  my_rxn->add_product ("N" ,chem_mixture.species_list_map().at(Antioch::Species::N) ,1);
  my_rxn->initialize();
  my_rxn->print();

//
  Scalar rate_reversible = my_rxn->compute_rate_of_progress(molar_densities,T,P0_RT,h_RT_minus_s_R);

//irreversible reaction
  my_rxn->set_reversibility(false);
  Scalar rate_irreversible = my_rxn->compute_rate_of_progress(molar_densities,T,P0_RT,h_RT_minus_s_R);

  int return_flag = 0;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;

  if( abs( (rate_reversible - net_rate)/net_rate) > tol )
    {
      std::cout << "Error: Mismatch in " << testname << " rate reversibility." << std::endl
                << std::setprecision(25)
                << "net_rate_antioch(T) = " << rate_reversible << std::endl
                << "net_rate_theory =    " << net_rate << std::endl
                << "relative error: " << abs( (rate_reversible - net_rate)/net_rate) << std::endl
                << "tolerance: " << tol << std::endl;

      return_flag = 1;
    }

  if( abs( (rate_irreversible - rate_fwd)/rate_fwd) > tol )
    {
      std::cout << "Error: Mismatch in " << testname << " rate reversibility." << std::endl
                << std::setprecision(25)
                << "irr_rate_antioch(T) = " << rate_reversible << std::endl
                << "irr_rate_theory =    " << rate_fwd << std::endl
                << "relative error: " << abs( (rate_irreversible - rate_fwd)/rate_fwd) << std::endl
                << "tolerance: " << tol << std::endl;

      return_flag = 1;
    }

  if( abs( (abs(rate_reversible - rate_irreversible) - rate_back)/rate_back) > tol )
    {
      std::cout << "Error: Mismatch in " << testname << " rate reversibility." << std::endl
                << std::setprecision(25)
                << "irrev_rate_antioch(T) - rev_rate_antioch(T) = " << rate_irreversible - rate_reversible << std::endl
                << "back_rate_theory =    " << rate_back << std::endl
                << "relative error: " << abs( ((rate_reversible - rate_irreversible) - rate_back)/rate_back) << std::endl
                << "tolerance: " << tol << std::endl;

      return_flag = 1;
    }

  delete my_rxn;

  return return_flag;
}

int main()
{
  return (tester<double>("double") ||
          tester<long double>("long double") ||
          tester<float>("float"));
}
