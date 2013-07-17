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
//
// $Id: arrhenius_rate_unit.C 38747 2013-04-17 23:26:39Z splessis $
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// C++
#include <limits>
// Antioch
#include "antioch/kooij_rate.h"

template <typename Scalar>
int tester()
{
  using std::abs;
  using std::exp;
  using std::pow;

  const Scalar Cf = 1.4;
  const Scalar eta = 1.2;
  const Scalar Ea = 5.0;

  Antioch::KooijRate<Scalar> kooij_rate(Cf,eta,Ea,1.,1.);

  int return_flag = 0;
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;

  for(Scalar T = 300.1; T <= 2500.1; T += 10.)
  {
  
    const Scalar rate_exact = Cf*pow(T,eta)*exp(-Ea/T);
    const Scalar derive_exact = exp(-Ea/T) * pow(T,eta) * Cf/T * (Ea/T + eta );

    Scalar rate1 = kooij_rate(T);
    Scalar deriveRate1 = kooij_rate.derivative(T);
    Scalar rate;
    Scalar deriveRate;

    kooij_rate.rate_and_derivative(T,rate,deriveRate);

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

int main()
{
  return (tester<double>() ||
          tester<long double>() ||
          tester<float>());
}
