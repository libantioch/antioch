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
#include "antioch/berthelot_rate.h"

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
int test_values(const Scalar & Cf, const Scalar & D, const Antioch::BerthelotRate<Scalar> & berthelot_rate)
{
  using std::abs;
  using std::exp;
  int return_flag = 0;

  for(Scalar T = 300.1; T <= 2500.1; T += 10.)
  {
    const Scalar rate_exact = Cf*exp(D*T);
    const Scalar derive_exact = D * Cf * exp(D*T);
    Antioch::KineticsConditions<Scalar> cond(T);

//KineticsConditions
    Scalar rate = berthelot_rate(cond);
    Scalar deriveRate = berthelot_rate.derivative(cond);

    return_flag = check_rate_and_derivative(rate_exact,derive_exact,rate,deriveRate,T) || return_flag;

    berthelot_rate.rate_and_derivative(cond,rate,deriveRate);

    return_flag = check_rate_and_derivative(rate_exact,derive_exact,rate,deriveRate,T) || return_flag;

// T
    rate = berthelot_rate(T);
    deriveRate = berthelot_rate.derivative(T);

    return_flag = check_rate_and_derivative(rate_exact,derive_exact,rate,deriveRate,T) || return_flag;

    berthelot_rate.rate_and_derivative(T,rate,deriveRate);

    return_flag = check_rate_and_derivative(rate_exact,derive_exact,rate,deriveRate,T) || return_flag;

  }
  return return_flag;
}

template <typename Scalar>
int tester()
{
  Scalar Cf = 1.4;
  Scalar D  = -5.0;

  Antioch::BerthelotRate<Scalar> berthelot_rate(Cf,D);

  int return_flag = test_values(Cf,D,berthelot_rate);

  Cf = 1e-7;
  D = 1e-3;
  berthelot_rate.set_Cf(Cf);
  berthelot_rate.set_D(D);
  return_flag = test_values(Cf,D,berthelot_rate) || return_flag;

  
  Cf = 2.1e-11;
  D = -1.2;
  std::vector<Scalar> values(2);
  values[0] = Cf;
  values[1] = D;
  berthelot_rate.reset_coefs(values);
  return_flag = test_values(Cf,D,berthelot_rate) || return_flag;

  return return_flag;
}

int main()
{
  return (tester<double>() ||
          tester<long double>() ||
          tester<float>());
}
