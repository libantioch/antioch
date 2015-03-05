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
// $Id: arrhenius_rate_unit.C 38747 2013-04-17 23:26:39Z splessis $
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// C++
#include <limits>
#include <vector>
// Antioch

#include "antioch/vector_utils_decl.h"

#include "antioch/kinetics_parsing.h"
//#include "antioch/physical_constants.h"
//#include "antioch/units.h"

#include "antioch/vector_utils.h"



template <typename Scalar>
int test_values(const Scalar & Cf, const Scalar & eta, const Scalar & Ea, const Scalar & D, const Scalar & Tref, const Scalar & R, const Antioch::KineticsType<Scalar> & rate_base)
{
  using std::abs;
  using std::exp;
  using std::pow;

  int return_flag = 0;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;

  for(Scalar T = 300.1L; T <= 2500.1L; T += 10.L)
  {
  
    const Scalar rate_exact = Cf*pow(T/Tref,eta)*exp(-Ea/(R*T) + D * T);
    const Scalar derive_exact = exp(-Ea/(R*T) + D * T) * pow(T/Tref,eta) * Cf * (Ea/(R*T*T) + eta/T + D );

    Scalar rate1 = rate_base(T);
    Scalar deriveRate1 = rate_base.derivative(T);
    Scalar rate;
    Scalar deriveRate;

    rate_base.compute_rate_and_derivative(T,rate,deriveRate);

    if( abs( (rate1 - rate_exact)/rate_exact ) > tol )
      {
          std::cout << std::scientific << std::setprecision(16)
                    << "Error: Mismatch in rate values." << std::endl
                    << "T = " << T << " K" << std::endl
                    << "rate(T) = " << rate1 << std::endl
                    << "rate_exact = " << rate_exact << std::endl
                    << "on rate " << rate_base << std::endl;

          return_flag = 1;
      }
    if( abs( (rate - rate_exact)/rate_exact ) > tol )
      {
          std::cout << std::scientific << std::setprecision(16)
                    << "Error: Mismatch in rate values." << std::endl
                    << "T = " << T << " K" << std::endl
                    << "rate(T) = " << rate << std::endl
                    << "rate_exact = " << rate_exact << std::endl
                    << "on rate " << rate_base << std::endl;

          return_flag = 1;
      }
    if( abs( (deriveRate1 - derive_exact)/derive_exact ) > tol )
      {
          std::cout << std::scientific << std::setprecision(16)
                    << "Error: Mismatch in rate derivative values." << std::endl
                    << "T = " << T << " K" << std::endl
                    << "drate_dT(T) = " << deriveRate1 << std::endl
                    << "derive_exact = " << derive_exact << std::endl
                    << "on rate " << rate_base << std::endl;

          return_flag = 1;
     }
    if( abs( (deriveRate - derive_exact)/derive_exact ) > tol )
      {
          std::cout << std::scientific << std::setprecision(16)
                    << "Error: Mismatch in rate derivative values." << std::endl
                    << "T = " << T << " K" << std::endl
                    << "drate_dT(T) = " << deriveRate << std::endl
                    << "derive_exact = " << derive_exact << std::endl
                    << "on rate " << rate_base << std::endl;

          return_flag = 1;
     }

     if(return_flag)break;
  }
  return return_flag;
}

template <typename Scalar>
int tester()
{

  const Scalar zero(0.L);
  Scalar Cf = 1.4L;
  Scalar eta = 1.2L;
  Scalar Ea = 298.0L;
  Scalar D = 0.05L;
  Scalar Tref = 1.L;
  Scalar R = 1.L;

/// building only here

   // constant
  std::vector<Scalar> coeffs(1,Cf); 

  Antioch::KineticsType<Scalar> * kin_base = Antioch::build_rate<Scalar>(coeffs,Antioch::KineticsModel::CONSTANT);

  int return_flag = test_values(Cf,zero,zero,zero,Tref,R,*kin_base);

   // Hercourt-Essen
  coeffs.resize(3);
  coeffs[0] = Cf;
  coeffs[1] = eta;
  coeffs[2] = Tref;

  delete kin_base;
  kin_base = Antioch::build_rate<Scalar>(coeffs,Antioch::KineticsModel::HERCOURT_ESSEN);
  
  return_flag = test_values(Cf,eta,zero,zero,Tref,R,*kin_base) || return_flag;

   // Arrhenius
  coeffs.resize(3);
  coeffs[0] = Cf;
  coeffs[1] = Ea;
  coeffs[2] = R;

  delete kin_base;
  kin_base = Antioch::build_rate<Scalar>(coeffs,Antioch::KineticsModel::ARRHENIUS);
  
  return_flag = test_values(Cf,zero,Ea,zero,Tref,R,*kin_base) || return_flag;

   // Berthelot
  coeffs.resize(2);
  coeffs[0] = Cf;
  coeffs[1] = D;

  delete kin_base;
  kin_base = Antioch::build_rate<Scalar>(coeffs,Antioch::KineticsModel::BERTHELOT);
  
  return_flag = test_values(Cf,zero,zero,D,Tref,R,*kin_base) || return_flag;

   // Berthelot Hercourt-Essen
  coeffs.resize(4);
  coeffs[0] = Cf;
  coeffs[1] = eta;
  coeffs[2] = D;
  coeffs[3] = Tref;

  delete kin_base;
  kin_base = Antioch::build_rate<Scalar>(coeffs,Antioch::KineticsModel::BHE);
  
  return_flag = test_values(Cf,eta,zero,D,Tref,R,*kin_base) || return_flag;

   // Kooij
  coeffs.resize(5);
  coeffs[0] = Cf;
  coeffs[1] = eta;
  coeffs[2] = Ea;
  coeffs[3] = Tref;
  coeffs[4] = R;

  delete kin_base;
  kin_base = Antioch::build_rate<Scalar>(coeffs,Antioch::KineticsModel::KOOIJ);
  
  return_flag = test_values(Cf,eta,Ea,zero,Tref,R,*kin_base) || return_flag;

   // Van't Hoff
  coeffs.resize(6);
  coeffs[0] = Cf;
  coeffs[1] = eta;
  coeffs[2] = Ea;
  coeffs[3] = D;
  coeffs[4] = Tref;
  coeffs[5] = R;

  delete kin_base;
  kin_base = Antioch::build_rate<Scalar>(coeffs,Antioch::KineticsModel::VANTHOFF);
  
  return_flag = test_values(Cf,eta,Ea,D,Tref,R,*kin_base) || return_flag;


//// reset
  Scalar Cf_reset = 1e-7L;
  Scalar eta_reset = 1.5L;
  Scalar Ea_reset = 36000.L;
  Scalar D_reset = -5.e-2L;
  Scalar Tref_reset = 298.;
  Scalar R_reset = Antioch::Constants::R_universal<Scalar>() * Antioch::Units<Scalar>("cal").get_SI_factor();

  coeffs[0] = Cf_reset;
  coeffs[1] = eta_reset;
  coeffs[2] = Ea_reset;
  coeffs[3] = D_reset;
  coeffs[4] = Tref_reset;
  coeffs[5] = R_reset;
  Antioch::reset_rate(*kin_base,coeffs);

  return_flag = test_values(Cf_reset,eta_reset,Ea_reset,D_reset,Tref_reset,R_reset,*kin_base) || return_flag;

  return return_flag;
}

int main()
{
  return (tester<double>() ||
          tester<long double>() ||
          tester<float>());
}
