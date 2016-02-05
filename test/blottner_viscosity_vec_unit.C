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



#include <valarray>

#ifdef ANTIOCH_HAVE_EIGEN
#include "Eigen/Dense"
#endif

#ifdef ANTIOCH_HAVE_METAPHYSICL
#include "metaphysicl/numberarray.h"
#endif

#ifdef ANTIOCH_HAVE_VEXCL
#include "vexcl/vexcl.hpp"
#endif

// Declare metaprogramming overloads before they're used
#include "antioch/eigen_utils_decl.h"
#include "antioch/metaphysicl_utils_decl.h"
#include "antioch/valarray_utils_decl.h"
#include "antioch/vexcl_utils_decl.h"

// C++
#include <cmath>
#include <iostream>
#include <limits>

// Antioch
#include "antioch/blottner_viscosity.h"

#include "antioch/eigen_utils.h"
#include "antioch/metaphysicl_utils.h"
#include "antioch/valarray_utils.h"
#include "antioch/vexcl_utils.h"

#ifdef ANTIOCH_HAVE_GRVY
#include "grvy.h"

GRVY::GRVY_Timer_Class gt;
#endif


template <typename Scalar, typename PairScalars>
int test_viscosity( const PairScalars mu, const PairScalars mu_exact,
		    const Scalar tol, const std::string& testname )
{
  using std::abs;

  int return_flag = 0;

  const PairScalars rel_error = abs( (mu - mu_exact)/mu_exact);

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

template <typename PairScalars>
int vectester(const PairScalars& example, const std::string& testname)
{
  typedef typename Antioch::value_type<PairScalars>::type Scalar;

  const Scalar a = 3.14e-3L;
  const Scalar b = 2.71e-2L;
  const Scalar c = 42.0e-5L;

  Antioch::BlottnerViscosity<Scalar> mu(a,b,c);

  std::cout << mu << std::endl;

  PairScalars T = example;
  PairScalars mu_exact = example;
  PairScalars mu_exact2 = example;

  for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
    {
      T[2*tuple]   = 1521.2L;
      T[2*tuple+1] = 1621.2L;

      // bc gives
      mu_exact[2*tuple]   = 0.1444222341677025337305172031086891L;
      mu_exact[2*tuple+1] = 0.1450979382180072302532592937776388L;

      mu_exact2[2*tuple]   = .1221724955488799960527696821225472L;
      mu_exact2[2*tuple+1] = .1224428450807678499433510473203746L;
    }

  int return_flag = 0;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 2;

#ifdef ANTIOCH_HAVE_GRVY
  gt.BeginTimer(testname);
#endif

  const PairScalars mu_T = mu(T);

#ifdef ANTIOCH_HAVE_GRVY
  gt.EndTimer(testname);
#endif

  return_flag = test_viscosity( mu_T, mu_exact, tol, testname );
  
  const Scalar a2 = 1e-3L;
  const Scalar b2 = 2e-2L;
  const Scalar c2 = 3e-5L;

  mu.reset_coeffs( a2, b2, c2 );

#ifdef ANTIOCH_HAVE_GRVY
  gt.BeginTimer(testname);
#endif

  const PairScalars mu2_T = mu(T);

#ifdef ANTIOCH_HAVE_GRVY
  gt.EndTimer(testname);
#endif

  return_flag = test_viscosity( mu2_T, mu_exact2, tol, testname );

  return return_flag;
}

int main()
{
  int returnval = 0;

  returnval = returnval ||
    vectester (std::valarray<float>(2*ANTIOCH_N_TUPLES), "valarray<float>");
  returnval = returnval ||
    vectester (std::valarray<double>(2*ANTIOCH_N_TUPLES), "valarray<double>");
  returnval = returnval ||
    vectester (std::valarray<long double>(2*ANTIOCH_N_TUPLES), "valarray<ld>");
#ifdef ANTIOCH_HAVE_EIGEN
  returnval = returnval ||
    vectester (Eigen::Array<float, 2*ANTIOCH_N_TUPLES, 1>(), "Eigen::ArrayXf");
  returnval = returnval ||
    vectester (Eigen::Array<double, 2*ANTIOCH_N_TUPLES, 1>(), "Eigen::ArrayXd");
  returnval = returnval ||
    vectester (Eigen::Array<long double, 2*ANTIOCH_N_TUPLES, 1>(), "Eigen::ArrayXld");
#endif
#ifdef ANTIOCH_HAVE_METAPHYSICL
  returnval = returnval ||
    vectester (MetaPhysicL::NumberArray<2*ANTIOCH_N_TUPLES, float> (0), "NumberArray<float>");
  returnval = returnval ||
    vectester (MetaPhysicL::NumberArray<2*ANTIOCH_N_TUPLES, double> (0), "NumberArray<double>");
  returnval = returnval ||
    vectester (MetaPhysicL::NumberArray<2*ANTIOCH_N_TUPLES, long double> (0), "NumberArray<ld>");
#endif
#ifdef ANTIOCH_HAVE_VEXCL
  vex::Context ctx_f (vex::Filter::All);
  if (!ctx_f.empty())
    returnval = returnval ||
      vectester (vex::vector<float> (ctx_f, 2*ANTIOCH_N_TUPLES), "vex::vector<float>");

  vex::Context ctx_d (vex::Filter::DoublePrecision);
  if (!ctx_d.empty())
    returnval = returnval ||
      vectester (vex::vector<double> (ctx_d, 2*ANTIOCH_N_TUPLES), "vex::vector<double>");
#endif

#ifdef ANTIOCH_HAVE_GRVY
  gt.Finalize();
  gt.Summarize();
#endif

  return returnval;
}
