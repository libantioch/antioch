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
#include <vector>
// Antioch
#include "antioch/vanthoff_rate.h"
#include "antioch/physical_constants.h"
#include "antioch/units.h"

template <typename Scalar>
int test_values(const Scalar & Cf, const Scalar & eta, const Scalar & Ea, const Scalar& D, const Scalar & Tref, const Scalar & R, const Antioch::VantHoffRate<Scalar> & vanthoff_rate)
{
  using std::abs;
  using std::exp;
  using std::pow;

  int return_flag = 0;
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;

  for(Scalar T = 300.1L; T <= 2500.1L; T += 10.L)
  {
    const Scalar rate_exact = Cf * pow(T/Tref,eta) * exp(-Ea/(R*T) + D*T);
    const Scalar derive_exact = Cf * pow(T/Tref,eta) * exp(-Ea/(R*T) + D*T) * (D + eta/T + Ea/(R*T*T));

    Scalar rate1 = vanthoff_rate(T);
    Scalar deriveRate1 = vanthoff_rate.derivative(T);
    Scalar rate;
    Scalar deriveRate;

    vanthoff_rate.rate_and_derivative(T,rate,deriveRate);

    if( abs( (rate1 - rate_exact)/rate_exact ) > tol )
      {
          std::cout << std::scientific << std::setprecision(16)
                    << "Error: Mismatch in rate values." << std::endl
                    << "T = " << T << " K" << std::endl
                    << "rate(T) = " << rate1 << std::endl
                    << "rate_exact = " << rate_exact << std::endl
                    << "Van't Hoff: " << vanthoff_rate << std::endl;

          return_flag = 1;
      }
    if( abs( (rate - rate_exact)/rate_exact ) > tol )
      {
          std::cout << std::scientific << std::setprecision(16)
                    << "Error: Mismatch in rate values." << std::endl
                    << "T = " << T << " K" << std::endl
                    << "rate(T) = " << rate << std::endl
                    << "rate_exact = " << rate_exact << std::endl
                    << "Van't Hoff: " << vanthoff_rate << std::endl;

          return_flag = 1;
      }
    if( abs( (deriveRate1 - derive_exact)/derive_exact ) > tol )
      {
          std::cout << std::scientific << std::setprecision(16)
                    << "Error: Mismatch in rate derivative values." << std::endl
                    << "T = " << T << " K" << std::endl
                    << "drate_dT(T) = " << deriveRate1 << std::endl
                    << "derive_exact = " << derive_exact << std::endl
                    << "Van't Hoff: " << vanthoff_rate << std::endl;

          return_flag = 1;
      }
    if( abs( (deriveRate - derive_exact)/derive_exact ) > tol )
      {
          std::cout << std::scientific << std::setprecision(16)
                    << "Error: Mismatch in rate derivative values." << std::endl
                    << "T = " << T << " K" << std::endl
                    << "drate_dT(T) = " << deriveRate << std::endl
                    << "derive_exact = " << derive_exact << std::endl
                    << "Van't Hoff: " << vanthoff_rate << std::endl;

          return_flag = 1;
      }
    if(return_flag)break;
  }

  return return_flag;
}

template <typename Scalar>
int tester()
{

  Scalar Cf = 1.4L;
  Scalar eta = 1.2L;
  Scalar Ea = 298.L;
  Scalar D = 2.50L;
  Scalar Tref = 1.L;
  Scalar R = 1.L;

  Antioch::VantHoffRate<Scalar> vanthoff_rate(Cf,eta,Ea,D,Tref,R);

  int return_flag = test_values(Cf,eta,Ea,D,Tref,R,vanthoff_rate);

  Cf = 1e-7L;
  eta = 0.8L;
  Ea = 36000.L;
  D = -5.0L;
  Tref = 298.;
  R = Antioch::Constants::R_universal<Scalar>() * Antioch::Units<Scalar>("cal").get_SI_factor();

  vanthoff_rate.set_Cf(Cf);
  vanthoff_rate.set_eta(eta);
  vanthoff_rate.set_Ea(Ea);
  vanthoff_rate.set_D(D);
  vanthoff_rate.set_Tref(Tref);
  vanthoff_rate.set_rscale(R);

  return_flag = test_values(Cf,eta,Ea,D,Tref,R,vanthoff_rate) || return_flag;

  Cf = 2.8e-7L;
  eta = 0.3L;
  Ea = 45000.L;
  D = -4.2L;
  std::vector<Scalar> values(4);
  values[0] = Cf;
  values[1] = eta;
  values[2] = Ea;
  values[3] = D;
  vanthoff_rate.reset_coefs(values);

  return_flag = test_values(Cf,eta,Ea,D,Tref,R,vanthoff_rate) || return_flag;

  Cf = 3.6e-11L;
  eta = 0.235L;
  Ea = 100000.L;
  D = 0.025L;
  Tref = 300.L;
  R = Antioch::Constants::R_universal<Scalar>();
  values.resize(6);
  values[0] = Cf;
  values[1] = eta;
  values[2] = Ea;
  values[3] = D;
  values[4] = Tref;
  values[5] = R;
  vanthoff_rate.reset_coefs(values);

  return_flag = test_values(Cf,eta,Ea,D,Tref,R,vanthoff_rate) || return_flag;


  return return_flag;
}

int main()
{
  return (tester<double>() ||
          tester<long double>() ||
          tester<float>());
}
