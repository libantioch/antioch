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
#include "antioch/constant_rate.h"

template <typename Scalar>
int check_rate_and_derivative(const Scalar & rate_exact, const Scalar & derive_exact,
                              const Scalar & rate, const Scalar & derive, const Scalar & T)
{
    const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 2;
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
int test_values(const Scalar & Cf, const Antioch::ConstantRate<Scalar> & constant_rate)
{
  int return_flag = 0;
  for(Scalar T = 300.1; T <= 2500.1; T += 10.)
  {

    const Scalar rate_exact = Cf;
    const Scalar derive_exact = 0.L;
    Antioch::KineticsConditions<Scalar> cond(T);

// KineticsConditions
    Scalar rate = constant_rate(cond);
    Scalar deriveRate = constant_rate.derivative(cond);

    return_flag = check_rate_and_derivative(rate_exact,derive_exact,rate,deriveRate,T) || return_flag;

    constant_rate.rate_and_derivative(cond,rate,deriveRate);

    return_flag = check_rate_and_derivative(rate_exact,derive_exact,rate,deriveRate,T) || return_flag;

// T
    rate = constant_rate(T);
    deriveRate = constant_rate.derivative(T);

    return_flag = check_rate_and_derivative(rate_exact,derive_exact,rate,deriveRate,T) || return_flag;

    constant_rate.rate_and_derivative(T,rate,deriveRate);

    return_flag = check_rate_and_derivative(rate_exact,derive_exact,rate,deriveRate,T) || return_flag;

  }
  return return_flag;
}

template <typename Scalar>
int tester()
{
  using std::abs;

  Scalar Cf = 1.4;

  Antioch::ConstantRate<Scalar> constant_rate(Cf);

  int return_flag = 0;

  return_flag = test_values(Cf,constant_rate) || return_flag;

  Cf = 1e-7;
  constant_rate.set_Cf(Cf);
  return_flag = test_values(Cf,constant_rate) || return_flag;

  Cf = 2.5e-11;
  std::vector<Scalar> values(1,Cf);
  constant_rate.reset_coefs(values);
  return_flag = test_values(Cf,constant_rate) || return_flag;

  return return_flag;
}

int main()
{
  return (tester<double>() ||
          tester<long double>() ||
          tester<float>());
}
