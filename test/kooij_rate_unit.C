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
#include "antioch/kooij_rate.h"
#include "antioch/physical_constants.h"
#include "antioch/units.h"

template <typename Scalar>
int check_rate_and_derivative(const Scalar & rate_exact, const Scalar & derive_exact, 
                              const Scalar & rate, const Scalar & derive, const Scalar & T)
{
    const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 2;
    int return_flag(0);
    if( abs( (rate - rate_exact)/rate_exact ) > tol )
      {
        std::cout << std::scientific << std::setprecision(16)
                  << "Error: Mismatch in rate values." << std::endl
                  << "T = " << T << " K" << std::endl
                  << "rate(T) = " << rate << std::endl
                  << "rate_exact = " << rate_exact << std::endl
                  << "relative difference = " <<  abs( (rate - rate_exact)/rate_exact ) << std::endl
                  << "tolerance = " <<  tol << std::endl;
 
        return_flag = 1;
      }
    if( abs( (derive - derive_exact)/derive_exact ) > tol )
      {
        std::cout << std::scientific << std::setprecision(16)
                  << "Error: Mismatch in rate derivative values." << std::endl
                  << "T = " << T << " K" << std::endl
                  << "drate_dT(T) = " << derive << std::endl
                  << "derive_exact = " << derive_exact << std::endl
                  << "relative difference = " <<  abs( (derive - derive_exact)/derive_exact ) << std::endl
                  << "tolerance = " <<  tol << std::endl;

        return_flag = 1;
      }

     return return_flag;
}

template <typename Scalar>
int test_values(const Scalar & Cf, const Scalar & eta, const Scalar & Ea, const Scalar & Tref, const Scalar & R, const Antioch::KooijRate<Scalar> & kooij_rate)
{
  using std::abs;
  using std::exp;
  using std::pow;

  int return_flag = 0;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;

  for(Scalar T = 300.1L; T <= 2500.1L; T += 10.L)
  {
  
    const Scalar rate_exact = Cf*pow(T/Tref,eta)*exp(-Ea/(R*T));
    const Scalar derive_exact = exp(-Ea/(R*T)) * pow(T/Tref,eta) * Cf * (Ea/(R*T*T) + eta/T );
    Antioch::KineticsConditions<Scalar> cond(T);

// KineticsConditions
    Scalar rate = kooij_rate(cond);
    Scalar deriveRate = kooij_rate.derivative(cond);

    return_flag = check_rate_and_derivative(rate_exact,derive_exact,rate,deriveRate,T) || return_flag;

    kooij_rate.rate_and_derivative(cond,rate,deriveRate);

    return_flag = check_rate_and_derivative(rate_exact,derive_exact,rate,deriveRate,T) || return_flag;

// T
    rate = kooij_rate(T);
    deriveRate = kooij_rate.derivative(T);

    return_flag = check_rate_and_derivative(rate_exact,derive_exact,rate,deriveRate,T) || return_flag;

    kooij_rate.rate_and_derivative(T,rate,deriveRate);

    return_flag = check_rate_and_derivative(rate_exact,derive_exact,rate,deriveRate,T) || return_flag;

  }
  return return_flag;
}

template <typename Scalar>
int tester()
{

  Scalar Cf = 1.4L;
  Scalar eta = 1.2L;
  Scalar Ea = 298.0L;
  Scalar Tref = 1.L;
  Scalar R = 1.L;

  Antioch::KooijRate<Scalar> kooij_rate(Cf,eta,Ea,Tref,R);

  int return_flag = test_values(Cf,eta,Ea,Tref,R,kooij_rate);

  Cf = 1e-7L;
  eta = 0.7L;
  Ea = 36000.0L;
  Tref = 298.;
  R = Antioch::Constants::R_universal<Scalar>() * Antioch::Units<Scalar>("cal").get_SI_factor();
  kooij_rate.set_Cf(Cf);
  kooij_rate.set_eta(eta);
  kooij_rate.set_Ea(Ea);
  kooij_rate.set_Tref(Tref);
  kooij_rate.set_rscale(R);

  return_flag = test_values(Cf,eta,Ea,Tref,R,kooij_rate) || return_flag;

  Cf = 2.1e-7L;
  eta = 1.1L;
  Ea = 23000.0L;
  std::vector<Scalar> values(3);
  values[0] = Cf;
  values[1] = eta;
  values[2] = Ea;
  kooij_rate.reset_coefs(values);

  return_flag = test_values(Cf,eta,Ea,Tref,R,kooij_rate) || return_flag;

  Cf = 2.1e-11L;
  eta = -0.01L;
  Ea = 100000.0L;
  R = Antioch::Constants::R_universal<Scalar>();
  Tref= 300.L;
  values.resize(5);
  values[0] = Cf;
  values[1] = eta;
  values[2] = Ea;
  values[3] = Tref;
  values[4] = R;
  kooij_rate.reset_coefs(values);

  return_flag = test_values(Cf,eta,Ea,Tref,R,kooij_rate) || return_flag;

  return return_flag;
}

int main()
{
  return (tester<double>() ||
          tester<long double>() ||
          tester<float>());
}
