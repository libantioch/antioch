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

#ifdef ANTIOCH_HAVE_VEXCL
#include "vexcl/vexcl.hpp"
#endif

// C++
#include <iostream>
#include <cmath>
#include <limits>

// Antioch

// Declare metaprogramming overloads before they're used
#include "antioch/eigen_utils_decl.h"
#include "antioch/metaphysicl_utils_decl.h"
#include "antioch/valarray_utils_decl.h"
#include "antioch/vector_utils_decl.h"
#include "antioch/vexcl_utils_decl.h"

#include "antioch/sutherland_viscosity.h"

#include "antioch/eigen_utils.h"
#include "antioch/metaphysicl_utils.h"
#include "antioch/valarray_utils.h"
#include "antioch/vexcl_utils.h"

#ifdef ANTIOCH_HAVE_GRVY
#include "grvy.h"

GRVY::GRVY_Timer_Class gt;
#endif

template <typename Scalar, typename PairScalars>
int test_viscosity( const PairScalars mu, const PairScalars mu_exact, const Scalar tol )
{
  using std::abs;

  int return_flag = 0;

  const PairScalars rel_error = abs( (mu - mu_exact)/mu_exact);

  if( Antioch::max(rel_error) > tol )
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


template <typename PairScalars>
int vectester(const PairScalars& example, const std::string& testname)
{
  typedef typename Antioch::value_type<PairScalars>::type Scalar;

  const Scalar mu_ref = 1.0e-3L;
  const Scalar T_ref = 300.0L;

  Antioch::SutherlandViscosity<Scalar> mu(mu_ref,T_ref);

  std::cout << mu << std::endl;

  PairScalars T = example;
  PairScalars mu_exact = example;
  PairScalars mu_exact2 = example;

  for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
    {
      T[2*tuple  ] = 1521.2L;
      T[2*tuple+1] = 1721.2L;

      // bc with scale=40 gives
      mu_exact[2*tuple  ] = .0325778060534850406481862157435995107036L;
      mu_exact[2*tuple+1] = .0353295183373055195000058747316029365368L;

      mu_exact2[2*tuple  ] = .0959985656417205050367745642313443587197L;
      mu_exact2[2*tuple+1] = .1047500160115581483776648561664869285592L;
    }

  int return_flag = 0;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 10;

#ifdef ANTIOCH_HAVE_GRVY
  gt.BeginTimer(testname);
#endif

  PairScalars muT = mu(T);

#ifdef ANTIOCH_HAVE_GRVY
  gt.EndTimer(testname);
#endif

  return_flag = test_viscosity( muT, mu_exact, tol );
  
  const Scalar mu_ref2 = 3.14159e-3L;
  const Scalar T_ref2 = 420.42L;

  mu.reset_coeffs(mu_ref2,T_ref2);

#ifdef ANTIOCH_HAVE_GRVY
  gt.BeginTimer(testname);
#endif

  muT = mu(T);

#ifdef ANTIOCH_HAVE_GRVY
  gt.EndTimer(testname);
#endif

  return_flag = test_viscosity( muT, mu_exact2, tol );

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
