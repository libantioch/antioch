//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

//Antioch
#include "antioch/vector_utils_decl.h"

#include "antioch/physical_constants.h"
#include "antioch/reaction_set.h"
#include "antioch/chemical_species.h"
#include "antioch/chemical_mixture.h"
#include "antioch/read_reaction_set_data.h"
#include "antioch/units.h"

#include "antioch/vector_utils.h"

//C++
#include <cmath>
#include <limits>

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

template<typename VectorScalar>
typename Antioch::value_type<VectorScalar>::type k_photo(const VectorScalar &solar_lambda, const VectorScalar &solar_irr, 
                                                         const VectorScalar &sigma_lambda, const VectorScalar &sigma_sigma)
{
  Antioch::SigmaBinConverter<VectorScalar> bin;
  VectorScalar sigma_rescaled;
  bin.y_on_custom_grid(sigma_lambda,sigma_sigma,solar_lambda,sigma_rescaled);

  typename Antioch::value_type<VectorScalar>::type _k(0.L);
  for(unsigned int il = 0; il < solar_irr.size() - 1; il++)
  {
     _k += sigma_rescaled[il] * solar_irr[il] * (solar_lambda[il+1] - solar_lambda[il]);
  }

  return _k;
}

template<typename Scalar>
int tester(const std::string &root_name)
{

  std::vector<std::string> species_str_list;
  species_str_list.push_back( "O2" );
  species_str_list.push_back( "OH" );
  species_str_list.push_back( "H2" );
  species_str_list.push_back( "H2O" );
  species_str_list.push_back( "H2O2" );
  species_str_list.push_back( "HO2" );
  species_str_list.push_back( "O" );
  species_str_list.push_back( "CH3" );
  species_str_list.push_back( "H" );
  unsigned int n_species = species_str_list.size();

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );
  Antioch::ReactionSet<Scalar> reaction_set( chem_mixture );
  Antioch::read_reaction_set_data_chemkin<Scalar>( root_name + "/test_parsing.chemkin", true, reaction_set );

//photochemistry set here
  std::vector<Scalar> hv,lambda;
  std::ifstream solar_flux(root_name + "/solar_flux.dat");
  std::string line;


//// the unit management here is tedious and useless, but it's got
//   all the steps, if ever someone needs a reference
  getline(solar_flux,line);
  Antioch::Units<Scalar> solar_wave("nm");
  Antioch::Units<Scalar> solar_irra("W/m2/nm");
  Antioch::Units<Scalar> i_unit = solar_irra - (Antioch::Constants::Planck_constant_unit<Scalar>() +  Antioch::Constants::light_celerity_unit<Scalar>() - solar_wave); //photons.s-1 = irradiance/(h*c/lambda)
  i_unit += Antioch::Units<Scalar>("nm"); //supress bin in unit calculations

  while(!solar_flux.eof())
  {
     Scalar l,i,di;
     solar_flux >> l >> i >> di;
     
     hv.push_back(i /(Antioch::Constants::Planck_constant<Scalar>() * Antioch::Constants::light_celerity<Scalar>() / l) // irr/(h*c/lambda): power -> number of photons.s-1
                                * i_unit.get_SI_factor()); //SI for cs, keep nm for bin
     lambda.push_back(l * solar_wave.factor_to_some_unit("nm")); //nm
     if(lambda.size() == 796)break;
  }
  solar_flux.close();

  std::vector<Scalar> CH4_s,CH4_lambda;
  std::ifstream CH4_file(root_name + "/CH4_hv_cs.dat");

  Scalar T = 2000.L;
  Scalar Tr = 1.;
  Antioch::Units<Scalar> unitA_m1("(cm3/mol)-1/s"),unitA_0("s-1"),unitA_1("cm3/mol/s"),unitA_2("(cm3/mol)2/s");

  Scalar Rcal = Antioch::Constants::R_universal<Scalar>() * Antioch::Constants::R_universal_unit<Scalar>().factor_to_some_unit("cal/mol/K");
  getline(CH4_file,line);

  Antioch::Units<Scalar> cs_input("cm2");
  Antioch::Units<Scalar> lambda_input("ang");
  Scalar factor_cs = cs_input.get_SI_factor() / lambda_input.factor_to_some_unit("nm");
  while(!CH4_file.eof())
  {
     Scalar l,s;
     CH4_file >> l >> s;
     CH4_s.push_back(s * factor_cs);
     CH4_lambda.push_back(l * lambda_input.factor_to_some_unit("nm"));
     if(CH4_s.size() == 137)break;
  }
  CH4_file.close();

  Antioch::ParticleFlux<std::vector<Scalar> > photons(lambda,hv);
  reaction_set.set_particle_flux(&photons);

//
  // Molar densities
  std::vector<Scalar> molar_densities(n_species,5e-4);
  Scalar tot_dens((Scalar)n_species * 5e-4);

///Elementary, + Kooij
  std::vector<Scalar> k;
  Scalar A,b,Ea;
/*
! Hessler, J. Phys. Chem. A, 102:4517 (1998)
H+O2=O+OH                 3.547e+15 -0.406  1.6599E+4
*/
 A  = 3.547e15L * unitA_m1.get_SI_factor();
 b  = -0.405L;
 Ea = 1.6599e4L; 
 k.push_back(Kooij(T,A,b,Ea,Tr,Rcal));

/*
! Sutherland et al., 21st Symposium, p. 929 (1986)
O+H2=H+OH                 0.508E+05  2.67  0.629E+04
*/
 k.push_back(Kooij(T,(Scalar)0.508e5,(Scalar)2.67,(Scalar)0.629e4,Tr,Rcal));

/*
! Michael and Sutherland, J. Phys. Chem. 92:3853 (1988)
H2+OH=H2O+H               0.216E+09  1.51  0.343E+04
*/
 k.push_back(Kooij(T,(Scalar)0.216e9,(Scalar)1.51,(Scalar)0.343e4,Tr,Rcal));

/*
! Sutherland et al., 23rd Symposium, p. 51 (1990)
O+H2O=OH+OH               2.97e+06   2.02  1.34e+4
*/
 k.push_back(Kooij(T,(Scalar)2.97e6,(Scalar)2.02,(Scalar)1.34e4,Tr,Rcal));

//! *************** H2-O2 Dissociation Reactions ******************
/*
! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986)
H2+M=H+H+M                4.577E+19 -1.40  1.0438E+05
   H2/2.5/ H2O/12/
*/
 k.push_back(Kooij(T,(Scalar)2.97e6,(Scalar)2.02,(Scalar)1.34e4,Tr,Rcal));

/*
! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986)
O+O+M=O2+M                6.165E+15 -0.50  0.000E+00
   H2/2.5/ H2O/12/
*/
 Scalar sum_eps = 5e-4L * (2.5L + 12.L);
 k.push_back(sum_eps * Kooij(T,(Scalar)6.165e15,(Scalar)-0.50,(Scalar)0.00,Tr,Rcal));

/*
! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986)
O+H+M=OH+M                4.714E+18 -1.00  0.000E+00
   H2/2.5/ H2O/12/
*/
 k.push_back(sum_eps * Kooij(T,(Scalar)4.714e18,(Scalar)-1.00,(Scalar)0.00,Tr,Rcal));

/*
! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986)
!H+OH+M=H2O+M              2.212E+22 -2.00  0.000E+00
H+OH+M=H2O+M               3.800E+22 -2.00  0.000E+00
   H2/2.5/ H2O/12/
*/
 k.push_back(sum_eps * Kooij(T,(Scalar)3.800e22,(Scalar)-2.00,(Scalar)0.00,Tr,Rcal));

/*
!************** Formation and Consumption of HO2******************

! Cobos et al., J. Phys. Chem. 89:342 (1985) for kinf
! Michael, et al., J. Phys. Chem. A, 106:5297 (2002) for k0

!******************************************************************************
*/
/*
! MAIN BATH GAS IS N2 (comment this reaction otherwise)
!
 H+O2(+M)=HO2(+M)      1.475E+12  0.60  0.00E+00
     LOW/6.366E+20  -1.72  5.248E+02/
     TROE/0.8  1E-30  1E+30/
     H2/2.0/ H2O/11./ O2/0.78/
*/
  Scalar k0   = Kooij(T,(Scalar)6.366e20,(Scalar)-1.72,(Scalar)5.248e2,Tr,Rcal);
  Scalar kinf = Kooij(T,(Scalar)1.475e12,(Scalar)0.60,(Scalar)0.00,Tr,Rcal);
  Scalar Pr = tot_dens * k0/kinf;
  Scalar Fc = FcentTroe(T,(Scalar)0.8,(Scalar)1e-30,(Scalar)1e30);
  k.push_back(k0 / (1./tot_dens + k0/kinf)  * FTroe(Fc,Pr));

/*
! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986) [modified]
HO2+H=H2+O2               1.66E+13   0.00   0.823E+03
*/
 k.push_back(Arrh(T,(Scalar)1.66e13,(Scalar)0.823e3,Rcal));

/*
! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986) [modified]
HO2+H=OH+OH               7.079E+13   0.00   2.95E+02
*/
 k.push_back(Arrh(T,(Scalar)7.079e13,(Scalar)2.95e2,Rcal));

/*
! Baulch et al., J. Phys. Chem. Ref Data, 21:411 (1992)
HO2+O=O2+OH               0.325E+14  0.00   0.00E+00
*/
 k.push_back(0.325e14);

/*
! Keyser, J. Phys. Chem. 92:1193 (1988)
HO2+OH=H2O+O2             2.890E+13  0.00 -4.970E+02
*/
 k.push_back(Arrh(T,(Scalar)2.890e13,(Scalar)-4.97e2,Rcal));

//! ***************Formation and Consumption of H2O2******************
/*
! Hippler et al., J. Chem. Phys. 93:1755 (1990)
HO2+HO2=H2O2+O2            4.200e+14  0.00  1.1982e+04
  DUPLICATE
HO2+HO2=H2O2+O2            1.300e+11  0.00 -1.6293e+3
  DUPLICATE
*/
 k.push_back(Arrh(T,(Scalar)4.200e14,(Scalar)1.1982e4,Rcal) + Arrh(T,(Scalar)1.300e11,(Scalar)-1.6293e3,Rcal));

/*
! Brouwer et al., J. Chem. Phys. 86:6171 (1987) for kinf
! Warnatz, J. in Combustion chemistry (1984) for k0
H2O2(+M)=OH+OH(+M)         2.951e+14   0.00  4.843E+04
  LOW/1.202E+17  0.00  4.55E+04/
  TROE/0.5 1E-30 1E+30/
  H2/2.5/ H2O/12/
*/
  k0   = Arrh(T,(Scalar)1.202e17,(Scalar)4.55e4,Rcal);
  kinf = Arrh(T,(Scalar)2.951e14,(Scalar)4.843e4,Rcal);
  Pr = tot_dens * k0/kinf;
  Fc = FcentTroe(T,(Scalar)0.5,(Scalar)1e-30,(Scalar)1e30);
  k.push_back(k0 / (1./tot_dens + k0/kinf)  * FTroe(Fc,Pr));
//

/*
! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986)
H2O2+H=H2O+OH             0.241E+14  0.00  0.397E+04
*/
 k.push_back(Arrh(T,(Scalar)0.241e14,(Scalar)0.397e4,Rcal));

/*
! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986)
H2O2+H=HO2+H2             0.482E+14  0.00  0.795E+04
*/
 k.push_back(Arrh(T,(Scalar)0.482e14,(Scalar)0.795e4,Rcal));

/*
! Tsang and Hampson, J. Phys. Chem. Ref. Data, 15:1087 (1986)
H2O2+O=OH+HO2             9.550E+06  2.00  3.970E+03
*/
 k.push_back(Kooij(T,(Scalar)9.550e6,(Scalar)2.00,(Scalar)3.970e3,Tr,Rcal));

/*
! Hippler and Troe, J. Chem. Phys. Lett. 192:333 (1992)
H2O2+OH=HO2+H2O           1.000E+12  0.00  0.000
    DUPLICATE
H2O2+OH=HO2+H2O           5.800E+14  0.00  9.557E+03
    DUPLICATE
*/
 k.push_back(Arrh(T,(Scalar)1.000e12,(Scalar)0.0) + Arrh(T,(Scalar)5.800e14,(Scalar)9.557e3,Rcal));

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;
  int return_flag(0);
  for(unsigned int ir = 0; ir < k.size(); ir++)
  {
     const Antioch::Reaction<Scalar> * reac = &reaction_set.reaction(ir);
     if(std::abs(k[ir] - reac->compute_forward_rate_coefficient(molar_densities,T))/k[ir] > tol)
     {
        std::cout << *reac << std::endl;
        std::cout << std::scientific << std::setprecision(16)
                  << "Error in kinetics comparison\n"
                  << "reaction #" << ir << "\n"
                  << "temperature: " << T << " K" << "\n"
                  << "theory: " << k[ir] << "\n"
                  << "calculated: " << reac->compute_forward_rate_coefficient(molar_densities,T) << "\n"
                  << "relative error = " << std::abs(k[ir] - reac->compute_forward_rate_coefficient(molar_densities,T))/k[ir] << "\n"
                  << "tolerance = " <<  tol
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
