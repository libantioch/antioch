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

// Antioch
#include "antioch/rotational_relaxation.h"
#include "antioch/metaprogramming.h"

// C++
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

const long double pi(3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825);

template <typename Scalar>
Scalar F(const Scalar & x)
{
  return 1.L + std::pow(pi,1.5)/2.L * std::sqrt(x) 
             + (pi*pi/4.L + 2.L) * x
             + std::pow(pi * x,1.5);
}

template <typename Scalar>
Scalar Z(const Scalar & T, const Scalar & eps_kb, const Scalar & z_298)
{
  return z_298 * F(eps_kb / 298.L) / F(eps_kb / T);
}


template <typename Scalar>
int tester()
{
  using std::abs;

  int return_flag = 0;
  const Scalar tol = (std::numeric_limits<Scalar>::epsilon() * 10 < 5e-17)?5e-17:
                      std::numeric_limits<Scalar>::epsilon() * 10;

  const Scalar eps_kb = 97.53L; // N2 value
  const Scalar z_298 = 4.0L; // N2 value

  Antioch::RotationalRelaxation<Scalar> rot(z_298,eps_kb);

  for(Scalar T = 300.1; T <= 2500.1; T += 10.)
  {
     Scalar z = rot(T);
     Scalar z_exact = Z(T,eps_kb,z_298);

    if( abs( (z - z_exact)/z_exact) > tol )
      {
          std::cout << std::scientific << std::setprecision(16)
                    << "Error: Mismatch in rotational relaxation values." << std::endl
                    << " T = " << T << std::endl
                    << " z = " << z << std::endl
                    << " z_exact = " << z_exact << std::endl
                    << " relative error = " << std::abs(z - z_exact)/z_exact << std::endl
                    << " tolerance = " << tol << std::endl;

          return_flag = 1;
     }
  }

  return return_flag;
}

int main()
{
  return (tester<double>() ||
          tester<long double>() ||
          tester<float>());
}
