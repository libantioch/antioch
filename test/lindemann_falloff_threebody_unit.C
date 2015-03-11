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

#include "antioch/kinetics_conditions.h"
#include "antioch/reaction.h"
#include "antioch/falloff_threebody_reaction.h"
#include "antioch/lindemann_falloff.h"

template <typename Scalar>
int tester(const std::string & type)
{
  using std::abs;
  using std::exp;
  using std::pow;


//values for 2 CH3 (+M) <=> C2H6 (+M) for the Kooij model, Ds are made up

  const Scalar Cf1 = 1.135e36L * 1e6L * 1e-12L; //(cm3/mol)^2/s -> kmol -> m3
  const Scalar beta1 = 1.246L; //true value is -5.246
  const Scalar Ea1 = 1704.8L / 1.9858775L; //cal/mol
  const Scalar D1 = -4e-2L; // K^-1
  const Scalar Cf2 = 6.22e16L * 1e3L * 1e-12L; //cm3/mol/s -> kmol -> m3
  const Scalar beta2 = -1.174L;
  const Scalar Ea2 = 635.8L / 1.9858775L; //cal/mol
  const Scalar D2 = -5e-3L;

  const std::string equation("A + B -> C + D");
  const unsigned int n_species(4);

  int return_flag = 0;
  std::vector<Scalar> mol_densities;
  mol_densities.push_back(1e-2);
  mol_densities.push_back(1e-2);
  mol_densities.push_back(1e-2);
  mol_densities.push_back(1e-2);

  std::vector<Scalar> epsilon;
  epsilon.push_back(2);
  epsilon.push_back(5);
  epsilon.push_back(8.5);
  epsilon.push_back(40);

  Scalar M = epsilon[0] * mol_densities[0];
  for(unsigned int i = 1; i < n_species; i++)
  {
     M += epsilon[i] * mol_densities[i];
  }

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 70;
  std::cout << type << ", tolerance = " << tol;
  Scalar max_diff(-1.L);

  for(Scalar T = 300.1; T <= 1500.1; T += 10.)
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
    rate_exact = k0 / (1.L/M + k0/kinf);
    derive_exact = rate_exact * (dk0_dT/k0 - dk0_dT/(kinf/M + k0) + k0 * dkinf_dT/(pow(kinf,2)/M + k0*kinf));
                  
    for(unsigned int i = 0; i < n_species; i++)
    {
      derive_dX_exact[i] = epsilon[i] * rate_exact/(M + pow(M,2)*k0/kinf);
    }

    Antioch::FalloffThreeBodyReaction<Scalar,Antioch::LindemannFalloff<Scalar> > * fall_reaction = 
        new Antioch::FalloffThreeBodyReaction<Scalar,Antioch::LindemannFalloff<Scalar> >(n_species,equation,true,Antioch::ReactionType::LINDEMANN_FALLOFF_THREE_BODY,kin_mod);

    fall_reaction->add_forward_rate(rate_kinetics1);
    fall_reaction->add_forward_rate(rate_kinetics2);
    for(unsigned int s = 0; s < n_species; s++)
    {
        fall_reaction->set_efficiency("",s,epsilon[s]);
    }

    Antioch::KineticsConditions<Scalar,std::vector<Scalar> > cond(T);
    Scalar rate1 = fall_reaction->compute_forward_rate_coefficient(mol_densities,cond);
    Scalar rate;
    Scalar drate_dT;
    std::vector<Scalar> drate_dx;
    drate_dx.resize(n_species);
    fall_reaction->compute_forward_rate_coefficient_and_derivatives(mol_densities,cond,rate,drate_dT,drate_dx);

    for(unsigned int i = 0; i < n_species; i++)
    {
      Scalar diff = abs( (drate_dx[i] - derive_dX_exact[i])/derive_dX_exact[i] );
      if(diff > max_diff)max_diff = diff;
      if( diff > tol )
      {
          std::cout << std::scientific << std::setprecision(16)
                    << "\nError: Mismatch in rate values." << std::endl
                    << "Kinetics model (see enum) " << kin_mod << std::endl
                    << "species " << i << std::endl
                    << "drate_dc(T) = " << drate_dx[i] << std::endl
                    << "drate_dc_exact = " << derive_dX_exact[i] << std::endl
                    << "relative error = " << abs((drate_dx[i] - derive_dX_exact[i])/derive_dX_exact[i]) << std::endl
                    << "tolerance = " << tol << std::endl;
          return_flag = 1;
      }
    }
    Scalar diff = abs( (rate1 - rate_exact)/rate_exact );
    if(diff > max_diff)max_diff = diff;
    if( diff > tol )
      {
        std::cout << std::scientific << std::setprecision(16)
                  << "\nError: Mismatch in rate values." << std::endl
                  << "Kinetics model (see enum) " << kin_mod << std::endl
                  << "T = " << T << " K" << std::endl
                  << "rate(T) = " << rate1 << std::endl
                  << "rate_exact = " << rate_exact << std::endl
                  << "relative error = " << abs((rate1 - rate_exact)/rate_exact) << std::endl
                  << "tolerance = " << tol << std::endl;

        return_flag = 1;
      }
    diff = abs( (rate - rate_exact)/rate_exact );
    if(diff > max_diff)max_diff = diff;
    if( diff > tol )
      {
        std::cout << std::scientific << std::setprecision(16)
                  << "\nError: Mismatch in rate values." << std::endl
                  << "Kinetics model (see enum) " << kin_mod << std::endl
                  << "T = " << T << " K" << std::endl
                  << "rate(T) = " << rate << std::endl
                  << "rate_exact = " << rate_exact << std::endl
                  << "relative error = " << abs((rate - rate_exact)/rate_exact) << std::endl
                  << "tolerance = " << tol << std::endl;

        return_flag = 1;
      }
    diff = abs( (drate_dT - derive_exact)/derive_exact );
    if(diff > max_diff)max_diff = diff;
    if( diff > tol )
      {
        std::cout << std::scientific << std::setprecision(16)
                  << "\nError: Mismatch in rate derivative values." << std::endl
                  << "Kinetics model (see enum) " << kin_mod << std::endl
                  << "T = " << T << " K" << std::endl
                  << "drate_dT(T) = " << drate_dT << std::endl
                  << "derive_exact = " << derive_exact << std::endl
                  << "relative error = " << abs((drate_dT - derive_exact)/derive_exact) << std::endl
                  << "tolerance = " << tol << std::endl;

        return_flag = 1;
      }

    delete fall_reaction;
    }
  }

  std::cout << " and maximum difference = " << max_diff << std::endl;
  return return_flag;
}

int main()
{
  return (tester<double>("double") ||
          tester<long double>("long double") ||
          tester<float>("float"));
}
