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
#include "antioch/falloff_reaction.h"

template <typename Scalar>
int tester()
{
  using std::abs;
  using std::exp;
  using std::pow;
  using std::log;

//values for 2 CH3 (+M) <=> C2H6 (+M) for the Kooij model, Ds are made up

  const Scalar Cf1 = 1.135e36 * 1e6 * 1e-12; //(cm3/mol)^2/s -> kmol -> m3
  const Scalar beta1 = 1.246; //true value is -5.246
  const Scalar Ea1 = 1704.8 / 1.9858775; //cal/mol
  const Scalar D1 = -4e-2; // K^-1
  const Scalar Cf2 = 6.22e16 * 1e3 * 1e-12; //cm3/mol/s -> kmol -> m3
  const Scalar beta2 = -1.174;
  const Scalar Ea2 = 635.8 / 1.9858775; //cal/mol
  const Scalar D2 = -5e-3;
  const Scalar alpha = 0.405;
  const Scalar T3 = 1120; //K
  const Scalar T1 = 69.6; //K
// T2 too high, negligible

  const std::string equation("A + B -> AB");
  const unsigned int n_species(3);

  int return_flag = 0;
  std::vector<Scalar> mol_densities;
  mol_densities.push_back(1e-2);
  mol_densities.push_back(1e-2);
  mol_densities.push_back(1e-2);
  Scalar M = mol_densities[0];
  for(unsigned int i = 1; i < n_species; i++)
  {
     M += mol_densities[i];
  }

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 1e6;

  for(Scalar T = 300.1; T <= 2500.1; T += 10.)
  {
    for(unsigned int ikinmod = 0; ikinmod < 6; ikinmod++)
    {

     Antioch::KineticsType<Scalar> *rate_kinetics1(NULL);
     Antioch::KineticsType<Scalar> *rate_kinetics2(NULL);
     Scalar k0,kinf,dk0_dT,dkinf_dT;
     Scalar rate_exact;
     Scalar derive_exact;
     std::vector<Scalar> derive_dX_exact;
     derive_dX_exact.resize(n_species);
//Troe pre-calculations
     Scalar Fcent = (1.L - alpha) * exp(-T/T3) + alpha * exp(-T/T1);
     Scalar dFcent_dT = (alpha - 1.L)/T3 * exp(-T/T3) - alpha/T1 * exp(-T/T1);
     Scalar dlog10Fcent_dT = Antioch::Constants::log10_to_log<Scalar>()/Fcent * dFcent_dT;
     Scalar n =   0.75L - 1.27L * Antioch::Constants::log10_to_log<Scalar>()*log(Fcent);
     Scalar c = - 0.40L - 0.67L * Antioch::Constants::log10_to_log<Scalar>()*log(Fcent);
     Scalar d = 0.14L;
     Scalar dn_dT = -1.27L * dlog10Fcent_dT;
     Scalar dc_dT = -0.67L * dlog10Fcent_dT;

     Antioch::KineticsModel::KineticsModel kin_mod;

    switch(ikinmod)
    {
      case 0:
      {
        kin_mod = Antioch::KineticsModel::HERCOURT_ESSEN;
        rate_kinetics1 =  new Antioch::HercourtEssenRate<Scalar>(Cf1,beta1,1.);
        rate_kinetics2 =  new Antioch::HercourtEssenRate<Scalar>(Cf2,beta2,1.);
        k0 = Cf1 * pow(T,beta1); kinf = Cf2 * pow(T,beta2);
        dk0_dT = Cf1 * pow (T,beta1) * beta1/T; dkinf_dT = Cf2 * pow (T,beta2) * beta2/T;
        break;
      }
      case 1:
      {
        kin_mod = Antioch::KineticsModel::BERTHELOT;
        rate_kinetics1 = new Antioch::BerthelotRate<Scalar>(Cf1,D1);
        rate_kinetics2 = new Antioch::BerthelotRate<Scalar>(Cf2,D2);
        k0 = Cf1 * exp(D1*T); kinf = Cf2 * exp(D2*T);
        dk0_dT = Cf1 * exp(D1*T) * D1; dkinf_dT = Cf2 * exp(D2*T) * D2;
        break;
      }
      case 2:
      {
        kin_mod = Antioch::KineticsModel::ARRHENIUS;
        rate_kinetics1 = new Antioch::ArrheniusRate<Scalar>(Cf1,Ea1,1.);
        rate_kinetics2 = new Antioch::ArrheniusRate<Scalar>(Cf2,Ea2,1.);
        k0 = Cf1 * exp(-Ea1/T); kinf = Cf2 * exp(-Ea2/T);
        dk0_dT = Cf1 * exp(-Ea1/T) * Ea1/pow(T,2); dkinf_dT = Cf2 * exp(-Ea2/T) * Ea2/pow(T,2);
        break;
      }
      case 3:
      {
        kin_mod = Antioch::KineticsModel::BHE;
        rate_kinetics1 = new Antioch::BerthelotHercourtEssenRate<Scalar>(Cf1,beta1,D1,1.);
        rate_kinetics2 = new Antioch::BerthelotHercourtEssenRate<Scalar>(Cf2,beta2,D2,1.);
        k0 = Cf1 * pow(T,beta1) * exp(D1 * T); kinf = Cf2 * pow(T,beta2) * exp(D2 * T);
        dk0_dT = Cf1 * pow(T,beta1) * exp(D1 * T) * (beta1/T + D1); dkinf_dT = Cf2 * pow(T,beta2) * exp(D2 * T) * (beta2/T + D2);
        break;
      }
      case 4:
      {
        kin_mod = Antioch::KineticsModel::KOOIJ;
        rate_kinetics1 = new Antioch::KooijRate<Scalar>(Cf1,beta1,Ea1,1.,1.);
        rate_kinetics2 = new Antioch::KooijRate<Scalar>(Cf2,beta2,Ea2,1.,1.);
        k0 = Cf1 * pow(T,beta1) * exp(-Ea1/T); kinf = Cf2 * pow(T,beta2) * exp(-Ea2/T);
        dk0_dT = Cf1 * pow(T,beta1) * exp(-Ea1/T) * (beta1/T + Ea1/pow(T,2)); dkinf_dT = Cf2 * pow(T,beta2) * exp(-Ea2/T) * (beta2/T + Ea2/pow(T,2));
        break;
      }
      case 5:
      {
        kin_mod = Antioch::KineticsModel::VANTHOFF;
        rate_kinetics1 = new Antioch::VantHoffRate<Scalar>(Cf1,beta1,Ea1,D1,1.,1.);
        rate_kinetics2 = new Antioch::VantHoffRate<Scalar>(Cf2,beta2,Ea2,D2,1.,1.);
        k0 = Cf1 * pow(T,beta1) * exp(-Ea1/T + D1 * T); kinf = Cf2 * pow(T,beta2) * exp(-Ea2/T + D2 * T);
        dk0_dT = Cf1 * pow(T,beta1) * exp(-Ea1/T + D1 * T) * (D1 + beta1/T + Ea1/pow(T,2));
        dkinf_dT = Cf2 * pow(T,beta2) * exp(-Ea2/T + D2 * T) * (D2 + beta2/T + Ea2/pow(T,2)); 
        break;
      }
    }
    // Troe calculations
    Scalar Pr = M * k0/kinf;
    Scalar dPr_dT = Pr * (dk0_dT/k0 - dkinf_dT/kinf);
    Scalar log10Pr = Antioch::Constants::log10_to_log<Scalar>() * log(Pr);
    Scalar dlog10Pr_dT = Antioch::Constants::log10_to_log<Scalar>() / Pr * dPr_dT;
    std::vector<Scalar> dlog10Pr_dX(n_species,0);
    for(unsigned int i = 0; i < n_species; i++)
    {
        dlog10Pr_dX[i] = Antioch::Constants::log10_to_log<Scalar>()/M;
    }
    Scalar logF = log(Fcent)/(1.L + pow(((log10Pr + c)/(n - d*(log10Pr + c) )),2));
    Scalar dlogF_dT = logF * (dlog10Fcent_dT / Fcent 
                                  - 2.L *pow((log10Pr + c)/(n - d * (log10Pr + c)),2) 
                                    * ((dlog10Pr_dT + dc_dT)/(log10Pr + c) -
                                       (dn_dT - d * (dlog10Pr_dT + dc_dT))/(n - d * (log10Pr + c))
                                      )
                                    / (1.L + pow((log10Pr + c)/(n - d * (log10Pr + c)),2)) 
                                 );
    Scalar F = exp(logF);
    Scalar dF_dT = F * dlogF_dT;
    std::vector<Scalar> dF_dX;
    dF_dX.resize(n_species);
    for(unsigned int i = 0; i < n_species; i++)
    {
        dF_dX[i] = - F * logF * logF/log(Fcent) * dlog10Pr_dX[i] * (1.L - 1.L/(n - d * (log10Pr + c))) * (log10Pr + c);
    }

    rate_exact = k0 / (1.L/M + k0/kinf);

    derive_exact = rate_exact * (dk0_dT/k0 - dk0_dT/(kinf/M + k0) + k0 * dkinf_dT/(kinf*(kinf/M + k0))) * F 
                 + dF_dT * rate_exact;

    for(unsigned int i = 0; i < n_species; i++)
    {
      derive_dX_exact[i] = rate_exact/(M + pow(M,2)*k0/kinf) * F + dF_dX[i] * rate_exact;
    }

    rate_exact *= F;

    Antioch::FalloffReaction<Scalar,Antioch::TroeFalloff<Scalar> > * fall_reaction = 
        new Antioch::FalloffReaction<Scalar,Antioch::TroeFalloff<Scalar> >(n_species,equation,true,Antioch::ReactionType::TROE_FALLOFF,kin_mod);
    fall_reaction->add_forward_rate(rate_kinetics1);
    fall_reaction->add_forward_rate(rate_kinetics2);
    fall_reaction->F().set_alpha(alpha);
    fall_reaction->F().set_T1(T1);
    fall_reaction->F().set_T3(T3);
    Scalar rate1 = fall_reaction->compute_forward_rate_coefficient(mol_densities,T);
    Scalar rate;
    Scalar drate_dT;
    std::vector<Scalar> drate_dx;
    drate_dx.resize(n_species);
    fall_reaction->compute_forward_rate_coefficient_and_derivatives(mol_densities,T,rate,drate_dT,drate_dx);

    for(unsigned int i = 0; i < n_species; i++)
    {
      if( abs( (drate_dx[i] - derive_dX_exact[i])/derive_dX_exact[i] ) > tol )
      {
          std::cout << std::scientific << std::setprecision(16)
                    << "Error: Mismatch in rate values." << std::endl
                    << "Kinetics model (see enum) " << kin_mod << std::endl
                    << "species " << i << std::endl
                    << "T = " << T << " K" << std::endl
                    << "drate_dc(T) = " << drate_dx[i] << std::endl
                    << "drate_dc_exact = " << derive_dX_exact[i] << std::endl
                    << "tolerance is " << tol << std::endl
                    << "criteria is " << abs( (drate_dx[i] - derive_dX_exact[i])/derive_dX_exact[i]) << std::endl
                    << "parameters are\nrate: " << rate_exact << "\n and  " << rate_exact << std::endl
                    << "total density: " << M << " species " << mol_densities[i] << std::endl;
          return_flag = 1;
      }
    }
    if(abs( (rate1 - rate_exact)/rate_exact ) > tol ) 
      {
        std::cout << std::scientific << std::setprecision(16)
                  << "Error: Mismatch in rate values." << std::endl
                  << "Kinetics model (see enum) " << kin_mod << std::endl
                  << "T = " << T << " K" << std::endl
                  << "rate(T) = " << rate1 << std::endl
                  << "rate_exact = " << rate_exact << std::endl
                  << "tolerance is " << tol << std::endl
                  << "criteria is " << abs( (rate1 - rate_exact)/rate_exact ) << std::endl;

        return_flag = 1;
      }
    if( abs( (rate - rate_exact)/rate_exact ) > tol )
      {
        std::cout << std::scientific << std::setprecision(16)
                  << "Error: Mismatch in rate values." << std::endl
                  << "Kinetics model (see enum) " << kin_mod << std::endl
                  << "T = " << T << " K" << std::endl
                  << "rate(T) = " << rate << std::endl
                  << "rate_exact = " << rate_exact << std::endl
                  << "tolerance is " << tol << std::endl
                  << "criteria is " << abs( (rate - rate_exact)/rate_exact ) << std::endl;

        return_flag = 1;
      }
    if( abs( (drate_dT - derive_exact)/derive_exact ) > tol )
      {
        std::cout << std::scientific << std::setprecision(16)
                  << "Error: Mismatch in rate derivative values." << std::endl
                  << "Kinetics model (see enum) " << kin_mod << std::endl
                  << "T = " << T << " K" << std::endl
                  << "drate_dT(T) = " << drate_dT << std::endl
                  << "derive_exact = " << derive_exact << std::endl
                  << "tolerance is " << tol << std::endl
                  << "criteria is " << abs( (drate_dT - derive_exact)/derive_exact ) << std::endl;

        return_flag = 1;
      }
    delete fall_reaction;
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
