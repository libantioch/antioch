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

// C++
#include <limits>
// Antioch
#include "antioch/vector_utils.h"

#include "antioch/reaction.h"
#include "antioch/elementary_reaction.h"

template <typename Scalar>
int tester()
{
  using std::abs;
  using std::exp;
  using std::pow;

  const Scalar Cf = 1.4;
  const Scalar Ea = 5.0;
  const Scalar beta = 1.2;
  const Scalar D = 2.5;

  const std::string equation("A + B -> C + D");
  const unsigned int n_species(4);

  int return_flag = 0;
  std::vector<Scalar> mol_densities;
  mol_densities.push_back(1e-2);
  mol_densities.push_back(1e-2);
  mol_densities.push_back(1e-2);
  mol_densities.push_back(1e-2);

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;

  for(Scalar T = 300.1; T <= 2500.1; T += 10.)
  {

    const Antioch::KineticsConditions<Scalar> conditions(T);

    for(unsigned int ikinmod = 0; ikinmod < 6; ikinmod++)
    {

     Antioch::KineticsType<Scalar> *rate_kinetics(NULL);
     Scalar rate_exact;
     Scalar derive_exact;
     Antioch::KineticsModel::KineticsModel kin_mod;

    switch(ikinmod)
    {
      case 0:
      {
        kin_mod = Antioch::KineticsModel::HERCOURT_ESSEN;
        rate_kinetics =  new Antioch::HercourtEssenRate<Scalar>(Cf,beta,1.);
        rate_exact = Cf * pow(T,beta);
        derive_exact = Cf * pow (T,beta) * beta/T;
        break;
      }
      case 1:
      {
        kin_mod = Antioch::KineticsModel::BERTHELOT;
        rate_kinetics = new Antioch::BerthelotRate<Scalar>(Cf,D);
        rate_exact = Cf * exp(D*T);
        derive_exact = Cf * exp(D*T) * D;
        break;
      }
      case 2:
      {
        kin_mod = Antioch::KineticsModel::ARRHENIUS;
        rate_kinetics = new Antioch::ArrheniusRate<Scalar>(Cf,Ea,1.);
        rate_exact = Cf * exp(-Ea/T);
        derive_exact = Cf * exp(-Ea/T) * Ea/pow(T,2);
        break;
      }
      case 3:
      {
        kin_mod = Antioch::KineticsModel::BHE;
        rate_kinetics = new Antioch::BerthelotHercourtEssenRate<Scalar>(Cf,beta,D,1.);
        rate_exact = Cf * pow(T,beta) * exp(D * T);
        derive_exact = Cf * pow(T,beta) * exp(D * T) * (beta/T + D);
        break;
      }
      case 4:
      {
        kin_mod = Antioch::KineticsModel::KOOIJ;
        rate_kinetics = new Antioch::KooijRate<Scalar>(Cf,beta,Ea,1.,1.);
        rate_exact = Cf * pow(T,beta) * exp(-Ea/T);
        derive_exact = Cf * pow(T,beta) * exp(-Ea/T) * (beta/T + Ea/pow(T,2));
        break;
      }
      case 5:
      {
        kin_mod = Antioch::KineticsModel::VANTHOFF;
        rate_kinetics = new Antioch::VantHoffRate<Scalar>(Cf,beta,Ea,D,1.,1.);
        rate_exact = Cf * pow(T,beta) * exp(-Ea/T + D * T);
        derive_exact =  Cf * pow(T,beta) * exp(-Ea/T + D * T) * (D + beta/T + Ea/pow(T,2));
        break;
      }
    }

    Antioch::ElementaryReaction<Scalar> * elem_reaction = new Antioch::ElementaryReaction<Scalar>(n_species,equation,true,kin_mod);
    elem_reaction->add_forward_rate(rate_kinetics);
    Scalar rate1 = elem_reaction->compute_forward_rate_coefficient(mol_densities,conditions);
    Scalar rate;
    Scalar drate_dT;
    std::vector<Scalar> drate_dx;
    drate_dx.resize(n_species);
    elem_reaction->compute_forward_rate_coefficient_and_derivatives(mol_densities,conditions,rate,drate_dT,drate_dx);

    if( abs( (rate1 - rate_exact)/rate_exact ) > tol )
      {
        std::cout << std::scientific << std::setprecision(16)
                  << "Error: Mismatch in rate values." << std::endl
                  << "Kinetics model (see enum) " << kin_mod << std::endl
                  << "T = " << T << " K" << std::endl
                  << "rate(T) = " << rate1 << std::endl
                  << "rate_exact = " << rate_exact << std::endl;

        return_flag = 1;
      }
    if( abs( (rate - rate_exact)/rate_exact ) > tol )
      {
        std::cout << std::scientific << std::setprecision(16)
                  << "Error: Mismatch in rate values." << std::endl
                  << "Kinetics model (see enum) " << kin_mod << std::endl
                  << "T = " << T << " K" << std::endl
                  << "rate(T) = " << rate << std::endl
                  << "rate_exact = " << rate_exact << std::endl;

        return_flag = 1;
      }
    if( abs( (drate_dT - derive_exact)/derive_exact ) > tol )
      {
        std::cout << std::scientific << std::setprecision(16)
                  << "Error: Mismatch in rate derivative values." << std::endl
                  << "Kinetics model (see enum) " << kin_mod << std::endl
                  << "T = " << T << " K" << std::endl
                  << "drate_dT(T) = " << drate_dT << std::endl
                  << "derive_exact = " << derive_exact << std::endl;

        return_flag = 1;
      }
    delete elem_reaction;
    if(return_flag)return return_flag;
    }
  }

  return return_flag;
}

int main()
{
  return (tester<double>() ||
          tester<long double>() ||
          tester<float>());
}
