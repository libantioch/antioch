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
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// C++
#include <iostream>
#include <cmath>

// Antioch
#include "antioch/blottner_viscosity.h"

template <typename Scalar>
int test_viscosity( const Scalar mu, const Scalar mu_exact, const Scalar tol )
{
  using std::abs;

  int return_flag = 0;

  const double rel_error = abs( (mu - mu_exact)/mu_exact);

  if( rel_error  > tol )
    {
      std::cerr << "Error: Mismatch in viscosity" << std::endl
		<< "mu(T)    = " << mu << std::endl
		<< "mu_exact = " << mu_exact << std::endl
		<< "rel_error = " << rel_error << std::endl
		<< "tol = " << tol << std::endl;
      return_flag = 1;
    }

  return return_flag;
}

template <typename Scalar>
int tester()
{
  const Scalar a = 3.14e-3;
  const Scalar b = 2.71e-2;
  const Scalar c = 42.0e-5;

  Antioch::BlottnerViscosity<Scalar> mu(a,b,c);

  std::cout << mu << std::endl;

  const Scalar T = 1521.2;

  // octave gives
  const Scalar mu_exact = 0.144422234167703;

  int return_flag = 0;

  const Scalar tol = 1.0e-14;

  return_flag = test_viscosity( mu(T), mu_exact, tol );
  
  const Scalar a2 = 1e-3;
  const Scalar b2 = 2e-2;
  const Scalar c2 = 3e-5;

  mu.reset_coeffs( a2, b2, c2 );

  // octave gives
  const Scalar mu_exact2 = 0.122172495548880;

  return_flag = test_viscosity( mu(T), mu_exact2, tol );

  return return_flag;
}

int main()
{
  return (tester<double>() ||
          tester<long double>() ||
          tester<float>());
}
