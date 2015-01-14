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
#include "antioch/threebody_reaction.h"

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
  std::vector<Scalar> species_eff;
  species_eff.push_back(0.4);
  species_eff.push_back(0.54);
  species_eff.push_back(0.5);
  species_eff.push_back(1.0);
  Scalar TB_term;
  Antioch::set_zero(TB_term);
  for(unsigned int i = 0; i < n_species; i++)
  {
     TB_term += species_eff[i] * mol_densities[i];
  }


  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;

  for(Scalar T = 300.1; T <= 2500.1; T += 10.)
  {

    const Antioch::KineticsConditions<Scalar> conditions(T);

    for(unsigned int ikinmod = 0; ikinmod < 6; ikinmod++)
    {

     Antioch::KineticsType<Scalar> *rate_kinetics(NULL);
     Scalar rate_exact;
     Scalar derive_exact;
     std::vector<Scalar> derive_dX_exact;
     derive_dX_exact.resize(n_species);
     Antioch::KineticsModel::KineticsModel kin_mod;

    switch(ikinmod)
    {
      case 0:
      {
        kin_mod = Antioch::KineticsModel::HERCOURT_ESSEN;
        rate_kinetics =  new Antioch::HercourtEssenRate<Scalar>(Cf,beta,1.);
        rate_exact = Cf * pow(T,beta) * TB_term;
        derive_exact = Cf * pow (T,beta) * beta/T * TB_term;
        for(unsigned int i = 0; i < n_species; i++)
        {
            derive_dX_exact[i] = Cf * pow(T,beta) * species_eff[i];
        }
        break;
      }
      case 1:
      {
        kin_mod = Antioch::KineticsModel::BERTHELOT;
        rate_kinetics = new Antioch::BerthelotRate<Scalar>(Cf,D);
        rate_exact = Cf * exp(D*T) * TB_term;
        derive_exact = Cf * exp(D*T) * D * TB_term;
        for(unsigned int i = 0; i < n_species; i++)
        {
            derive_dX_exact[i] = Cf * exp(D*T) * species_eff[i];
        }
        break;
      }
      case 2:
      {
        kin_mod = Antioch::KineticsModel::ARRHENIUS;
        rate_kinetics = new Antioch::ArrheniusRate<Scalar>(Cf,Ea,1.);
        rate_exact = Cf * exp(-Ea/T) * TB_term;
        derive_exact = Cf * exp(-Ea/T) * Ea/pow(T,2) * TB_term;
        for(unsigned int i = 0; i < n_species; i++)
        {
            derive_dX_exact[i] = Cf * exp(-Ea/T) * species_eff[i];
        }
        break;
      }
      case 3:
      {
        kin_mod = Antioch::KineticsModel::BHE;
        rate_kinetics = new Antioch::BerthelotHercourtEssenRate<Scalar>(Cf,beta,D,1.);
        rate_exact = Cf * pow(T,beta) * exp(D * T) * TB_term;
        derive_exact = Cf * pow(T,beta) * exp(D * T) * (beta/T + D) * TB_term;
        for(unsigned int i = 0; i < n_species; i++)
        {
            derive_dX_exact[i] = Cf * pow(T,beta) * exp(D * T) * species_eff[i];
        }
        break;
      }
      case 4:
      {
        kin_mod = Antioch::KineticsModel::KOOIJ;
        rate_kinetics = new Antioch::KooijRate<Scalar>(Cf,beta,Ea,1.,1.);
        rate_exact = Cf * pow(T,beta) * exp(-Ea/T) * TB_term;
        derive_exact = Cf * pow(T,beta) * exp(-Ea/T) * (beta/T + Ea/pow(T,2)) * TB_term;
        for(unsigned int i = 0; i < n_species; i++)
        {
            derive_dX_exact[i] = Cf * pow(T,beta) * exp(-Ea/T) * species_eff[i];
        }
        break;
      }
      case 5:
      {
        kin_mod = Antioch::KineticsModel::VANTHOFF;
        rate_kinetics = new Antioch::VantHoffRate<Scalar>(Cf,beta,Ea,D,1.,1.);
        rate_exact = Cf * pow(T,beta) * exp(-Ea/T + D * T) * TB_term;
        derive_exact =  Cf * pow(T,beta) * exp(-Ea/T + D * T) * (D + beta/T + Ea/pow(T,2)) * TB_term;
        for(unsigned int i = 0; i < n_species; i++)
        {
            derive_dX_exact[i] = Cf * pow(T,beta) * exp(-Ea/T + D * T) * species_eff[i];
        }
        break;
      }
    }

    Antioch::ThreeBodyReaction<Scalar> * TB_reaction = new Antioch::ThreeBodyReaction<Scalar>(n_species,equation,true,kin_mod);
    TB_reaction->add_forward_rate(rate_kinetics);
    for(unsigned int i = 0; i < n_species; i++)
    {
        TB_reaction->set_efficiency("",i,species_eff[i]);
    }
    Scalar rate1 = TB_reaction->compute_forward_rate_coefficient(mol_densities,conditions);
    Scalar rate;
    Scalar drate_dT;
    std::vector<Scalar> drate_dx;
    drate_dx.resize(n_species);
    TB_reaction->compute_forward_rate_coefficient_and_derivatives(mol_densities,conditions,rate,drate_dT,drate_dx);


    for(unsigned int i = 0; i < n_species; i++)
      {
      if( abs( (drate_dx[i] - derive_dX_exact[i])/derive_dX_exact[i] ) > tol )
        {
          std::cout << std::scientific << std::setprecision(16)
                    << "Error: Mismatch in derivative values." << std::endl
                    << "Kinetics model (see enum) " << kin_mod << std::endl
                    << "species " << i << std::endl
                    << "species molar densities " << mol_densities[i] << std::endl
                    << "species efficiency " << species_eff[i] << std::endl
                    << "drate_dX = " << drate_dx[i] << std::endl
                    << "drate_dX_exact = " << derive_dX_exact[i] << std::endl;
        return_flag = 1;
        }
      }
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
    delete TB_reaction;
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
