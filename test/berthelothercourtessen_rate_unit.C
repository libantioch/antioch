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
#include "antioch/berthelothercourtessen_rate.h"


template <typename Scalar>
int test_values(const Scalar & Cf, const Scalar & eta, const Scalar & D, const Scalar & Tref, const Antioch::BerthelotHercourtEssenRate<Scalar> & berthelothercourtessen_rate)
{
  using std::abs;
  using std::exp;
  using std::pow;
  int return_flag = 0;
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;

  for(Scalar T = 300.1L; T <= 2500.1L; T += 10.L)
  {
    const Scalar rate_exact = Cf*pow(T/Tref,eta)*exp(D*T);
    const Scalar derive_exact = Cf * pow(T/Tref,eta) * exp(D*T) * (D + eta/T);

    Scalar rate1 = berthelothercourtessen_rate(T);
    Scalar deriveRate1 = berthelothercourtessen_rate.derivative(T);
    Scalar rate;
    Scalar deriveRate;

    berthelothercourtessen_rate.rate_and_derivative(T,rate,deriveRate);

    if( abs( (rate1 - rate_exact)/rate_exact ) > tol )
      {
          std::cout << std::scientific << std::setprecision(16)
                    << "Error: Mismatch in rate values." << std::endl
                    << "T = " << T << " K" << std::endl
                    << "rate(T) = " << rate1 << std::endl
                    << "rate_exact = " << rate_exact << std::endl;

          return_flag = 1;
      }
    if( abs( (rate - rate_exact)/rate_exact ) > tol )
      {
          std::cout << std::scientific << std::setprecision(16)
                    << "Error: Mismatch in rate values." << std::endl
                    << "T = " << T << " K" << std::endl
                    << "rate(T) = " << rate << std::endl
                    << "rate_exact = " << rate_exact << std::endl;

          return_flag = 1;
      }
    if( abs( (deriveRate1 - derive_exact)/derive_exact ) > tol )
      {
          std::cout << std::scientific << std::setprecision(16)
                    << "Error: Mismatch in rate derivative values." << std::endl
                    << "T = " << T << " K" << std::endl
                    << "drate_dT(T) = " << deriveRate1 << std::endl
                   << "derive_exact = " << derive_exact << std::endl;

          return_flag = 1;
      }
    if( abs( (deriveRate - derive_exact)/derive_exact ) > tol )
      {
          std::cout << std::scientific << std::setprecision(16)
                    << "Error: Mismatch in rate derivative values." << std::endl
                    << "T = " << T << " K" << std::endl
                    << "drate_dT(T) = " << deriveRate << std::endl
                   << "derive_exact = " << derive_exact << std::endl;

          return_flag = 1;
      }

    if(return_flag)break;
  }
  return return_flag;

}


template <typename Scalar>
int tester()
{
  Scalar Cf  = 1.4L;
  Scalar eta = 1.2L;
  Scalar D   = -5.0L;
  Scalar Tref   = 1.0L;

  Antioch::BerthelotHercourtEssenRate<Scalar> berthelothercourtessen_rate(Cf,eta,D);

  int return_flag = test_values(Cf,eta,D,Tref,berthelothercourtessen_rate);

  Cf  = 1e-7L;
  eta = 0.6L;
  D   = 1e-3L;
  Tref   = 298.0L;
  berthelothercourtessen_rate.set_Cf(Cf);
  berthelothercourtessen_rate.set_eta(eta);
  berthelothercourtessen_rate.set_D(D);
  berthelothercourtessen_rate.set_Tref(Tref);

  return_flag = test_values(Cf,eta,D,Tref,berthelothercourtessen_rate) || return_flag;

  Cf  = 2.1e-7L;
  eta = 0.35L;
  D   = 8.1e-3L;
  std::vector<Scalar> values(3);
  values[0] = Cf;
  values[1] = eta;
  values[2] = D;
  berthelothercourtessen_rate.reset_coefs(values);

  return_flag = test_values(Cf,eta,D,Tref,berthelothercourtessen_rate) || return_flag;

  Cf  = 2.1e-11L;
  eta = -0.35L;
  D   = -2.1;
  Tref   = 300.L;
  values.resize(4);
  values[0] = Cf;
  values[1] = eta;
  values[2] = D;
  values[3] = Tref;
  berthelothercourtessen_rate.reset_coefs(values);

  return_flag = test_values(Cf,eta,D,Tref,berthelothercourtessen_rate) || return_flag;


  return return_flag;
}

int main()
{
  return (tester<double>() ||
          tester<long double>() ||
          tester<float>());
}
