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

#include "antioch_config.h"

#include <valarray>

#ifdef ANTIOCH_HAVE_EIGEN
#include "Eigen/Dense"
#endif

#ifdef ANTIOCH_HAVE_METAPHYSICL
#include "metaphysicl/numberarray.h"
#endif

// Declare metaprogramming overloads before they're used
#include "antioch/eigen_utils_decl.h"
#include "antioch/metaphysicl_utils_decl.h"
#include "antioch/valarray_utils_decl.h"

// C++
#include <cmath>
#include <iostream>
#include <limits>

// Antioch
#include "antioch/blottner_viscosity.h"

#include "antioch/eigen_utils.h"
#include "antioch/metaphysicl_utils.h"
#include "antioch/valarray_utils.h"

template <typename Scalar, typename PairScalars>
int test_viscosity( const PairScalars mu, const PairScalars mu_exact, const Scalar tol )
{
  int return_flag = 0;

  const PairScalars rel_error = std::abs( (mu - mu_exact)/mu_exact);

  if( Antioch::max(rel_error) > tol )
    {
      std::cerr << "Error: Mismatch in viscosity" << std::endl
		<< "mu(T)    = (" << mu[0] << ',' << mu[1] << ')' << std::endl
		<< "mu_exact = (" << mu_exact[0] << ',' << mu_exact[1] << ')' << std::endl
		<< "rel_error = (" << rel_error[0] << ',' << rel_error[1] << ')' << std::endl
		<< "tol = " << tol << std::endl;
      return_flag = 1;
    }

  return return_flag;
}

template <typename Scalar, typename PairScalars>
int vectester(const PairScalars& example)
{
  const Scalar a = 3.14e-3;
  const Scalar b = 2.71e-2;
  const Scalar c = 42.0e-5;

  Antioch::BlottnerViscosity<Scalar> mu(a,b,c);

  std::cout << mu << std::endl;

  PairScalars T = example;
  T[0] = 1521.2;
  T[1] = 1621.2;

  // bc gives
  PairScalars mu_exact = example;
  mu_exact[0] = 0.1444222341677025337305172031086891L;
  mu_exact[1] = 0.1450979382180072302532592937776388L;

  int return_flag = 0;

  // How are we getting such high error in the long double case?
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 600;

  return_flag = test_viscosity( mu(T), mu_exact, tol );
  
  const Scalar a2 = 1e-3;
  const Scalar b2 = 2e-2;
  const Scalar c2 = 3e-5;

  mu.reset_coeffs( a2, b2, c2 );

  // octave gives
  PairScalars mu_exact2 = example;
  mu_exact2[0] = .1221724955488799960527696821225472L;
  mu_exact2[1] = .1224428450807678499433510473203746L;

  return_flag = test_viscosity( mu(T), mu_exact2, tol );

  return return_flag;
}

int main()
{
  int returnval = 0;

  returnval = returnval ||
    vectester<float, std::valarray<float> >
      (std::valarray<float>(2));
  returnval = returnval ||
    vectester<double, std::valarray<double> >
      (std::valarray<double>(2));
  returnval = returnval ||
    vectester<long double, std::valarray<long double> >
      (std::valarray<long double>(2));
#ifdef ANTIOCH_HAVE_EIGEN
  returnval = returnval ||
    vectester<float, Eigen::Array2f>
      (Eigen::Array2f());
  returnval = returnval ||
    vectester<double, Eigen::Array2d>
      (Eigen::Array2d());
  returnval = returnval ||
    vectester<long double, Eigen::Array<long double, 2, 1> >
      (Eigen::Array<long double, 2, 1>());
#endif
#ifdef ANTIOCH_HAVE_METAPHYSICL
  returnval = returnval ||
    vectester<float, MetaPhysicL::NumberArray<2, float> > (0);
  returnval = returnval ||
    vectester<double, MetaPhysicL::NumberArray<2, double> > (0);
  returnval = returnval ||
    vectester<long double, MetaPhysicL::NumberArray<2, long double> > (0);
#endif

  return returnval;
}
