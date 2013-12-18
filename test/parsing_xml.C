//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
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

//Antioch
#include "antioch/physical_constants.h"
#include "antioch/reaction_set.h"
#include "antioch/chemical_species.h"
#include "antioch/chemical_mixture.h"
#include "antioch/read_reaction_set_data_xml.h"
#include "antioch/units.h"

//C++
#include <cmath>
#include <limits>

template<typename Scalar>
Scalar HE(const Scalar &T, const Scalar &Cf, const Scalar &eta, const Scalar &Tf = 1.)
{
  return Cf * std::pow(T/Tf,eta);
}
template<typename Scalar>
Scalar Bert(const Scalar &T, const Scalar &Cf, const Scalar &D)
{
  return Cf * std::exp(D * T);
}

template<typename Scalar>
Scalar BHE(const Scalar &T, const Scalar &Cf, const Scalar &eta, const Scalar &D, const Scalar &Tf = 1.)
{
  return Cf * std::pow(T/Tf,eta) * std::exp(D * T);
}

template<typename Scalar>
Scalar Arrh(const Scalar &T, const Scalar &Cf, const Scalar &Ea, const Scalar &R = Antioch::Constants::R_universal<Scalar>())
{
  return Cf * std::exp(-Ea /(R * T));
}

template<typename Scalar>
Scalar Kooij(const Scalar &T, const Scalar &Cf, const Scalar &eta, const Scalar &Ea, const Scalar &Tf = 1., 
             const Scalar &R = Antioch::Constants::R_universal<Scalar>())
{
  return Cf * std::pow(T/Tf,eta) * std::exp(-Ea /(R * T));
}

template<typename Scalar>
Scalar VH(const Scalar &T, const Scalar &Cf, const Scalar &eta, const Scalar &Ea, const Scalar &D, 
          const Scalar &Tf = 1., const Scalar &R = Antioch::Constants::R_universal<Scalar>())
{
  return Cf * std::pow(T/Tf,eta) * std::exp(-Ea /(R * T) + D * T);
}

template<typename Scalar>
Scalar FcentTroe(const Scalar &T, const Scalar &alpha, const Scalar &T3, const Scalar &T1, const Scalar &T2 = -1.)
{
  return (T2 < Scalar(0.))?(Scalar(1.) - alpha) * std::exp(-T/T3) + alpha * std::exp(-T/T1):
                           (Scalar(1.) - alpha) * std::exp(-T/T3) + alpha * std::exp(-T/T1) + std::exp(-T2/T);
}

template<typename Scalar>
Scalar coeffTroe(const Scalar &coef1, const Scalar &coef2, const Scalar &Fcent)
{
  return coef1 + coef2 * std::log10(Fcent);
}

template<typename Scalar>
Scalar cTroe(const Scalar &Fcent)
{
  Scalar c1(-0.40);
  Scalar c2(-0.67);
  return coeffTroe(c1,c2,Fcent);
}

template<typename Scalar>
Scalar nTroe(const Scalar &Fcent)
{
  Scalar c1(0.75);
  Scalar c2(-1.27);
  return coeffTroe(c1,c2,Fcent);
}

template<typename Scalar>
Scalar FTroe(const Scalar &Fcent, const Scalar &Pr)
{
  Scalar d(0.14L);
  return std::pow(10.L,std::log10(Fcent)/(Scalar(1.) + 
               std::pow( (std::log10(Pr) + cTroe(Fcent)) /
                                ( nTroe(Fcent) - d * (std::log10(Pr) + cTroe(Fcent)) ),2)));
}

template<typename Scalar>
int tester(const std::string &input_name)
{

  std::vector<std::string> species_str_list;
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "O2" );
  species_str_list.push_back( "N" );
  species_str_list.push_back( "O" );
  species_str_list.push_back( "NO" );
  species_str_list.push_back( "C" );
  species_str_list.push_back( "C2" );
  species_str_list.push_back( "CN" );
  unsigned int n_species = species_str_list.size();

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );
  Antioch::ReactionSet<Scalar> reaction_set( chem_mixture );
  Antioch::read_reaction_set_data_xml<Scalar>( input_name, true, reaction_set );


  Scalar T = 2000.L;
  Scalar Tr = 1.;
  Antioch::Units<Scalar> unitA_m1("(m3/kmol)-1/s"),unitA_0("s-1"),unitA_1("m3/kmol/s"),unitA_2("(m3/kmol)2/s");

  Scalar Rcal = Antioch::Constants::R_universal<Scalar>() * Antioch::Constants::R_universal_unit<Scalar>().factor_to_some_unit("cal/mol/K");
  // Molar densities
  std::vector<Scalar> molar_densities(n_species,5e-4);

///Elementary, + Kooij - Arrhenius conversion tested
  std::vector<Scalar> k;
  Scalar A,beta,Ea,D;

// N2 -> 2 N
  A    = 7e18 * unitA_0.get_SI_factor();
  beta = -1.6;
  k.push_back(HE(T,A,beta));

// O2 -> 2 O
  A = 2e18 * unitA_0.get_SI_factor();
  D = -5e-3;
  k.push_back(Bert(T,A,D));

//NO -> N + O
  A  = 5e12 * unitA_0.get_SI_factor();
  Ea = 149943.0;
  k.push_back(Arrh(T,A,Ea,Rcal));
  beta = 0.42;
  k.push_back(Kooij(T,A,beta,Ea,Tr,Rcal));

//N2 + O -> NO + N
  A = 5.7e9 * unitA_1.get_SI_factor();
  beta = 0.42;
  k.push_back(BHE(T,A,beta,D));

//NO + O -> NO + N
  A = 8.4e9 * unitA_1.get_SI_factor();
  beta = 0.4;
  Ea = 38526.0;
  k.push_back(Kooij(T,A,beta,Ea,Tr,Rcal));
  k.push_back(Arrh(T,A,Ea,Rcal));

//C2 -> 2 C
  A = 3.7e11 * unitA_0.get_SI_factor();
  beta = -0.42;
  Ea = 138812.8;
  k.push_back(Kooij(T,A,beta,Ea,Tr,Rcal));

//CN -> C + N
  A = 2.5e11 * unitA_0.get_SI_factor();
  beta = 0.40;
  Ea = 174240.9;
  D = 0.05;
  k.push_back(VH(T,A,beta,Ea,D,Tr,Rcal));

///Duplicate
  Scalar A2,beta2,Ea2,D2,A3,beta3,Ea3,D3,A4,Ea4;

// N2 -> 2 N
  A     = 7e18 * unitA_0.get_SI_factor();
  beta  = -1.6;
  A2    = 5e17 * unitA_0.get_SI_factor();
  beta2 = 0.5;
  A3    = 3e18 * unitA_0.get_SI_factor();
  beta3 = -0.6;
  k.push_back(HE(T,A,beta) + HE(T,A2,beta2) + HE(T,A3,beta3));

// O2 -> 2 O
  A  = 2e18 * unitA_0.get_SI_factor();
  D  = -5e-2;
  A2 = 2e+16 * unitA_0.get_SI_factor();
  D2 = 0.003;
  k.push_back(Bert(T,A,D) + Bert(T,A2,D2));

// NO -> N + O
  A   = 5e+12 * unitA_0.get_SI_factor();
  Ea  = 149943.0;
  A2  = 3.5e+10 * unitA_0.get_SI_factor();
  Ea2 = 1943.0;
  A3  = 1.5e+8 * unitA_0.get_SI_factor();
  Ea3 = 149.0;
  A4  = 5.5e+8 * unitA_0.get_SI_factor();
  Ea4 = 943.0;
  k.push_back(Arrh(T,A,Ea,Rcal) + Arrh(T,A2,Ea2,Rcal) + Arrh(T,A3,Ea3,Rcal) + Arrh(T,A4,Ea4,Rcal));

// N2 + O -> NO + N
  A     = 5.7e+9 * unitA_1.get_SI_factor();
  beta  = 0.42;
  D     = -5e-3;
  A2    = 7e+7 * unitA_1.get_SI_factor();
  beta2 = 0.5;
  D2    = 2.5e-5;
  k.push_back(BHE(T,A,beta,D) + BHE(T,A2,beta2,D2));

//NO + O -> NO + N
  A     = 8.4e+09 * unitA_1.get_SI_factor();
  beta  = 0.40;
  Ea    = 38526.0;
  A2    = 4e+07 * unitA_1.get_SI_factor();
  beta2 = 0.50;
  Ea2   = 40500.0;
  A3    = 5e+10 * unitA_1.get_SI_factor();
  beta3 = 0.10;
  Ea3   = 15000.0;
  k.push_back(Kooij(T,A,beta,Ea,Tr,Rcal) + Kooij(T,A2,beta2,Ea2,Tr,Rcal) + Kooij(T,A3,beta3,Ea3,Tr,Rcal));

//C2 -> 2 C
  A     = 3.7e+11 * unitA_0.get_SI_factor();
  beta  = -0.42;
  Ea    = 138812.8;
  A2    = 5.0e+10 * unitA_0.get_SI_factor();
  beta2 = 1.32;
  Ea2   = 150500.8;
  k.push_back(Kooij(T,A,beta,Ea,Tr,Rcal) + Kooij(T,A2,beta2,Ea2,Tr,Rcal));

//CN -> C + N
  A     = 2.5e+11 * unitA_0.get_SI_factor();
  beta  = 0.40;
  D     = -5e-3;
  Ea    = 174240.9;
  A2    = 5e+10 * unitA_0.get_SI_factor();
  beta2 = 0.50;
  D2    = -1.5e-2;
  Ea2   = 4240.9;
  A3    = 3.2e+10 * unitA_0.get_SI_factor();
  beta3 = 1.20;
  D3    = -2.5e-5;
  Ea3   = 174.9;
  k.push_back(VH(T,A,beta,Ea,D,Tr,Rcal) + VH(T,A2,beta2,Ea2,D2,Tr,Rcal) + VH(T,A3,beta3,Ea3,D3,Tr,Rcal));

//three body
// N2 -> 2 N
  A    = 7e18 * unitA_1.get_SI_factor();
  beta = -1.6;
  Ea   = 149943.0;
  k.push_back(HE(T,A,beta) * (Scalar(n_species) - 2. + 4.2857 + 4.2857) * 5e-4);

// O2 -> 2 O
  A = 2e18 * unitA_1.get_SI_factor();
  D = -5e-3;
  k.push_back(Bert(T,A,D) * (Scalar(n_species) - 2. + 5.0 + 5.0) * 5e-4);

//NO -> N + O
  A = 5e12 * unitA_1.get_SI_factor();
  k.push_back(Arrh(T,A,Ea,Rcal) * (Scalar(n_species) - 3. + 22.0 + 22.0 + 22.0) * 5e-4);

//N2 + O -> NO + N
  A = 5.7e9 * unitA_2.get_SI_factor();
  beta = 0.42;
  D = -5e-3;
  k.push_back(BHE(T,A,beta,D) * (Scalar(n_species) - 3. + 22.0 + 22.0 + 22.0) * 5e-4);

//NO + O -> NO + N
  A = 8.4e9 * unitA_2.get_SI_factor();
  beta = 0.4;
  Ea = 38526.0;
  k.push_back(Kooij(T,A,beta,Ea,Tr,Rcal) * (Scalar(n_species) - 3. + 22.0 + 22.0 + 22.0) * 5e-4);

//C2 -> 2 C
  A = 3.7e11 * unitA_1.get_SI_factor();
  beta = -0.42;
  Ea = 138812.8;
  k.push_back(Kooij(T,A,beta,Ea,Tr,Rcal) * Scalar(n_species) * 5e-4);

//CN -> C + N
  A = 2.5e11 * unitA_1.get_SI_factor();
  beta = 0.40;
  Ea = 729372.4;
  D = 5e-3;
  k.push_back(VH(T,A,beta,Ea,D,Tr,Rcal)  * Scalar(n_species) * 5e-4);
///Lindemann Falloff
// falloff is k(T,[M]) = k0*[M]/(1 + [M]*k0/kinf) * F = k0 * ([M]^-1 + k0 * kinf^-1)^-1 * F    
// F = 1

// N2 -> 2 N
  A     = 7e18 * unitA_1.get_SI_factor();
  beta  = -1.6;
  A2    = 5e15 * unitA_0.get_SI_factor();
  beta2 = 0.5;
  k.push_back(HE(T,A,beta) / (1./4e-3 + HE(T,A,beta)/HE(T,A2,beta2)) );

// O2 -> 2 O
  A  = 5e17 * unitA_1.get_SI_factor();
  D  = -2.5e-5;
  A2 = 2e18 * unitA_0.get_SI_factor();
  D2 = -5e-3;
  k.push_back(Bert(T,A,D) / (1./4e-3 + Bert(T,A,D)/Bert(T,A2,D2)) );

//NO -> N + O
  A   = 5.e+12 * unitA_1.get_SI_factor();
  Ea  = 149943.0;
  A2  = 3e+15 * unitA_0.get_SI_factor();
  Ea2 =  200000.0;
  k.push_back(Arrh(T,A,Ea,Rcal) / (1./4e-3 + Arrh(T,A,Ea,Rcal)/Arrh(T,A2,Ea2,Rcal)) );

//N2 + O -> NO + N
  A     = 5e+9 * unitA_2.get_SI_factor();
  beta  = 0.6;
  D     = -5e-4;
  A2    = 5.7e+9 * unitA_1.get_SI_factor();
  beta2 = -0.42;
  D2    = -5e-3;
  k.push_back(BHE(T,A,beta,D) / (1./4e-3 + BHE(T,A,beta,D)/BHE(T,A2,beta2,D2)) );

//NO + O -> NO + N
  A     = 8.4e+09 * unitA_2.get_SI_factor();
  beta  = 0.40;
  Ea    = 38526.0;
  A2    = 8.4e+05 * unitA_1.get_SI_factor();
  beta2 = 0.02;
  Ea2   = 3526.0;
  k.push_back(Kooij(T,A,beta,Ea,Tr,Rcal) / (1./4e-3 + Kooij(T,A,beta,Ea,Tr,Rcal)/Kooij(T,A2,beta2,Ea2,Tr,Rcal)) );

//C2 -> 2 C
  A     = 3.7e+11 * unitA_1.get_SI_factor();
  beta  = -0.42;
  Ea    = 138812.8;
  A2    = 3.7e+12 * unitA_0.get_SI_factor();
  beta2 = -0.52;
  Ea2   = 135000.8;
  k.push_back(Kooij(T,A,beta,Ea,Tr,Rcal) / (1./4e-3 + Kooij(T,A,beta,Ea,Tr,Rcal)/Kooij(T,A2,beta2,Ea2,Tr,Rcal)) );

//CN -> C + N
  A     = 5e+10 * unitA_1.get_SI_factor();
  beta  = -0.10;
  D     = 1.5e-3;
  Ea    = 150240.9;
  A2    = 2.5e+11 * unitA_0.get_SI_factor();
  beta2 = 0.40;
  D2    = -0.005;
  Ea2   = 174240.9;
  k.push_back(VH(T,A,beta,Ea,D,Tr,Rcal) / (1./4e-3 + VH(T,A,beta,Ea,D,Tr,Rcal)/VH(T,A2,beta2,Ea2,D2,Tr,Rcal)) );
//Troe falloff
//falloff is k(T,[M]) = k0*[M]/(1 + [M]*k0/kinf) * F = k0 * ([M]^-1 + k0 * kinf^-1)^-1 * F    
// F is complicated...
  Scalar Pr,k0,kinf;
  Scalar Fc,alpha,T1,T2,T3;
  alpha = 0.562;
  T1    = 5836;
  T2    = 8552;
  T3    = 91;
  Fc = FcentTroe(T,alpha,T3,T1,T2);

// N2 -> 2 N
  A     = 7.e+18 * unitA_1.get_SI_factor();
  beta  = -1.6;
  A2    = 5.e+15 * unitA_0.get_SI_factor();
  beta2 = 0.5;
  k0   = HE(T,A,beta);
  kinf = HE(T,A2,beta2);
  Pr = 4e-3 * k0/kinf;
  k.push_back(k0 / (1./4e-3 + k0/kinf)  * FTroe(Fc,Pr));

// O2 -> 2 O
  A  = 5e17 * unitA_1.get_SI_factor();
  D  = -2.5e-5;
  A2 = 2e18 * unitA_0.get_SI_factor();
  D2 = -5e-3;
  k0   = Bert(T,A,D);
  kinf = Bert(T,A2,D2);
  Pr = 4e-3 * k0/kinf;
  k.push_back(k0 / (1./4e-3 + k0/kinf)  * FTroe(Fc,Pr));

//NO -> N + O
  A   = 5.e+12 * unitA_1.get_SI_factor();
  Ea  = 149943.0;
  A2  = 3e+15 * unitA_0.get_SI_factor();
  Ea2 =  200000.0;
  k0    = Arrh(T,A,Ea,Rcal);
  kinf  = Arrh(T,A2,Ea2,Rcal);
  Pr = 4e-3 * k0/kinf;
  k.push_back(k0 / (1./4e-3 + k0/kinf)  * FTroe(Fc,Pr));

//N2 + O -> NO + N
  A     = 5e+9 * unitA_2.get_SI_factor();
  beta  = 0.6;
  D     = -5e-4;
  A2    = 5.7e+9 * unitA_1.get_SI_factor();
  beta2 = -0.42;
  D2    = -5e-3;
  k0    = BHE(T,A,beta,D); 
  kinf  = BHE(T,A2,beta2,D2);
  Pr = 4e-3 * k0/kinf;
  k.push_back(k0 / (1./4e-3 + k0/kinf)  * FTroe(Fc,Pr));

//NO + O -> NO + N
  A     = 8.4e+09 * unitA_2.get_SI_factor();
  beta  = 0.40;
  Ea    = 38526.0;
  A2    = 8.4e+05 * unitA_1.get_SI_factor();
  beta2 = 0.02;
  Ea2   = 3526.0;
  k0    = Kooij(T,A,beta,Ea,Tr,Rcal);
  kinf  = Kooij(T,A2,beta2,Ea2,Tr,Rcal);
  Pr = 4e-3 * k0/kinf;
  k.push_back(k0 / (1./4e-3 + k0/kinf)  * FTroe(Fc,Pr));

//C2 -> 2 C
  A     = 3.7e+11 * unitA_1.get_SI_factor();
  beta  = -0.42;
  Ea    = 138812.8;
  A2    = 3.7e+12 * unitA_0.get_SI_factor();
  beta2 = -0.52;
  Ea2   = 135000.8;
  k0    = Kooij(T,A,beta,Ea,Tr,Rcal); 
  kinf  = Kooij(T,A2,beta2,Ea2,Tr,Rcal);
  Pr = 4e-3 * k0/kinf;
  k.push_back(k0 / (1./4e-3 + k0/kinf)  * FTroe(Fc,Pr));

//CN -> C + N
  A     = 5e+10 * unitA_1.get_SI_factor();
  beta  = -0.10;
  D     = 1.5e-3;
  Ea    = 150240.9;
  A2    = 2.5e+11 * unitA_0.get_SI_factor();
  beta2 = 0.40;
  D2    = -0.005;
  Ea2   = 174240.9;
  k0    = VH(T,A,beta,Ea,D,Tr,Rcal); 
  kinf  = VH(T,A2,beta2,Ea2,D2,Tr,Rcal);
  Pr = 4e-3 * k0/kinf;
  k.push_back(k0 / (1./4e-3 + k0/kinf)  * FTroe(Fc,Pr));

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;
  int return_flag(0);
  for(unsigned int ir = 0; ir < k.size(); ir++)
  {
     const Antioch::Reaction<Scalar> * reac = &reaction_set.reaction(ir);
     if(std::abs(k[ir] - reac->compute_forward_rate_coefficient(molar_densities,T))/k[ir] > tol)
     {
        std::cout << std::scientific << std::setprecision(16)
                  << "Error in kinetics comparison\n"
                  << "reaction #: " << ir << "\n"
                  << "theory: " << k[ir] << "\n"
                  << "calculated: " << reac->compute_forward_rate_coefficient(molar_densities,T) << "\n"
                  << "relative error = " << std::abs(k[ir] - reac->compute_forward_rate_coefficient(molar_densities,T))/k[ir] << "\n"
                  <<  "tolerance = " <<  tol
                  << std::endl;
        return_flag = 1;
     }
  }

  return return_flag;
}

int main(int argc, char* argv[])
{
  return (tester<float>(std::string(argv[1])) ||
          tester<double>(std::string(argv[1])));/* ||
          tester<long double>(std::string(argv[1])));*/
}
