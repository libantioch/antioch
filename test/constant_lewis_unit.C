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
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// C++
#include <iostream>
#include <cmath>
#include <limits>

// Antioch
#include "antioch/constant_lewis_diffusivity.h"

template <typename Scalar>
int test_val( const Scalar val, const Scalar val_exact, const Scalar tol, const std::string& val_name )
{
  using std::abs;

  int return_flag = 0;

  const Scalar rel_error = abs( (val - val_exact)/val_exact);

  if( rel_error  > tol )
    {
      std::cerr << "Error: Mismatch in " << val_name << std::endl
                << val_name << "    = " << val << std::endl
                << val_name+"_exact = " << val_exact << std::endl
                << "rel_error = " << rel_error << std::endl
                << "tol = " << tol << std::endl;
      return_flag = 1;
    }

  return return_flag;
}

template <typename Scalar>
int tester()
{
  const Scalar Le = 1.4;

  Antioch::ConstantLewisDiffusivity<Scalar> diff( Le );

  const Scalar rho = 3.14;

  const Scalar cp = 2.71;

  const Scalar k = 42.0;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 2;

  const Scalar D = diff.D( rho, cp, k );

  const Scalar D_exact = k/(Le*rho*cp);

  return test_val( D, D_exact, tol, std::string("D") );
}

int main()
{
  return ( tester<double>()  ||
           tester<long double>() ||
           tester<float>()
           );
}
