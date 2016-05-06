//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014-2016 Paul T. Bauman, Benjamin S. Kirk,
//                         Sylvain Plessis, Roy H. Stonger
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

// C++
#include <limits>
#include <vector>
// Antioch
#include "antioch/hercourtessen_rate.h"

template <typename Scalar>
int check_rate_and_derivative(const Scalar & rate_exact, const Scalar & derive_exact,
                              const Scalar & rate, const Scalar & derive, const Scalar & T)
{
    const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;
    int return_flag(0);

    using std::abs;

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
int test_values(const Scalar & Cf, const Scalar & eta, const Scalar & Tref, const Antioch::HercourtEssenRate<Scalar> & hercourtessen_rate)
{
  using std::abs;
  using std::pow;
  int return_flag = 0;

  for(Scalar T = 300.1L; T <= 2500.1L; T += 10.L)
  {
    const Scalar rate_exact = Cf * pow(T/Tref,eta);
    const Scalar derive_exact = Cf * eta * pow(T/Tref,eta)/T;
    Antioch::KineticsConditions<Scalar> cond(T);


// KineticsConditions
    Scalar rate = hercourtessen_rate(cond);
    Scalar deriveRate = hercourtessen_rate.derivative(cond);

    return_flag = check_rate_and_derivative(rate_exact,derive_exact,rate,deriveRate,T) || return_flag;

    hercourtessen_rate.rate_and_derivative(cond,rate,deriveRate);

    return_flag = check_rate_and_derivative(rate_exact,derive_exact,rate,deriveRate,T) || return_flag;

// T
    rate = hercourtessen_rate(T);
    deriveRate = hercourtessen_rate.derivative(T);

    return_flag = check_rate_and_derivative(rate_exact,derive_exact,rate,deriveRate,T) || return_flag;

    hercourtessen_rate.rate_and_derivative(T,rate,deriveRate);

    return_flag = check_rate_and_derivative(rate_exact,derive_exact,rate,deriveRate,T) || return_flag;

  }
  return return_flag;
}

template <typename Scalar>
int tester()
{

  Scalar Cf = 1.4L;
  Scalar eta = 1.2L;
  Scalar Tref = 1.L;

  Antioch::HercourtEssenRate<Scalar> hercourtessen_rate(Cf,eta); // default

  int return_flag = test_values(Cf,eta,Tref,hercourtessen_rate);
  Cf = 1e-7L;
  eta = 0.7L;
  Tref = 300.L;
  hercourtessen_rate.set_Cf(Cf);
  hercourtessen_rate.set_eta(eta);
  hercourtessen_rate.set_Tref(Tref);
  return_flag = test_values(Cf,eta,Tref,hercourtessen_rate) || return_flag;

  Cf = 2.1e-11L;
  eta = -2.3L;
  std::vector<Scalar> values(2);
  values[0] = Cf;
  values[1] = eta;
  hercourtessen_rate.reset_coefs(values);
  return_flag = test_values(Cf,eta,Tref,hercourtessen_rate) || return_flag;

  Tref = 298.L;
  values.resize(3);
  values[0] = Cf;
  values[1] = eta;
  values[2] = Tref;
  hercourtessen_rate.reset_coefs(values);
  return_flag = test_values(Cf,eta,Tref,hercourtessen_rate) || return_flag;

  return return_flag;
}

int main()
{
  return (tester<double>() ||
          tester<long double>() ||
          tester<float>());
}
