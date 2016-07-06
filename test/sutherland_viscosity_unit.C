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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// C++
#include <iostream>
#include <cmath>
#include <limits>

// Antioch
#include "antioch/sutherland_viscosity.h"


template <typename Scalar>
int test_viscosity( const Scalar mu, const Scalar mu_exact, const Scalar tol )
{
  using std::abs;

  int return_flag = 0;

  const Scalar rel_error = abs( (mu - mu_exact)/mu_exact);

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
  const Scalar mu_ref = 1.0e-3L;
  const Scalar T_ref = 300.0L;

  Antioch::SutherlandViscosity<Scalar> mu(mu_ref,T_ref);

  std::cout << mu << std::endl;

  const Scalar T = 1521.2L;

  // bc with scale=40 gives
  const Scalar mu_exact = .0325778060534850406481862157435995107036L;

  int return_flag = 0;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 10;

  return_flag = test_viscosity( mu(T), mu_exact, tol );
  
  const Scalar mu_ref2 = 3.14159e-3L;
  const Scalar T_ref2 = 420.42L;

  mu.reset_coeffs(mu_ref2,T_ref2);

  // bc with scale=40 gives
  const Scalar mu_exact2 = .0959985656417205050367745642313443587197L;

  return_flag = test_viscosity( mu(T), mu_exact2, tol );

  return return_flag;
}


int main()
{
  return (tester<double>() ||
          tester<long double>() ||
          tester<float>());
}
