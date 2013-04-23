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
// $Id: kooij_rate_vec_unit.C 38747 2013-04-17 23:26:39Z splessis $
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// valarray has to be declared before Antioch or gcc can't find the
// right versions of exp() and pow() to use??

#include "antioch_config.h"

#include <valarray>

#ifdef ANTIOCH_HAVE_EIGEN
#include "Eigen/Dense"
#endif

#ifdef ANTIOCH_HAVE_METAPHYSICL
#include "metaphysicl/numberarray.h"
#endif

#include "antioch/kooij_rate.h"

#include "antioch/eigen_utils.h"
#include "antioch/metaphysicl_utils.h"
#include "antioch/valarray_utils.h"

#include <cmath>
#include <limits>

template <typename PairScalars>
int vectester(const PairScalars& example)
{
  using std::abs;
  using std::exp;
  using std::pow;

  typedef typename Antioch::value_type<PairScalars>::type Scalar;

  const Scalar Cf = 1.4;
  const Scalar eta = 1.2;
  const Scalar Ea = 5.0;

  Antioch::KooijRate<Scalar> kooij_rate(Cf,eta,Ea);

  // Construct from example to avoid resizing issues
  PairScalars T = example;
  T[0] = 1500.1;
  T[1] = 1600.1;
  
  const Scalar rate_exact0 = Cf*pow(Scalar(1500.1),eta)*exp(-Ea/1500.1);
  const Scalar rate_exact1 = Cf*pow(Scalar(1600.1),eta)*exp(-Ea/1600.1);

  int return_flag = 0;

  const PairScalars rate = kooij_rate(T);

  const Scalar tol = std::numeric_limits<Scalar>::epsilon()*10;

  if( abs( (rate[0] - rate_exact0)/rate_exact0 ) > tol )
    {
      std::cout << "Error: Mismatch in rate values." << std::endl
		<< "rate(T0)   = " << rate[0] << std::endl
		<< "rate_exact = " << rate_exact0 << std::endl
		<< "difference = " << rate[0] - rate_exact0 << std::endl;

      return_flag = 1;
    }

  if( abs( (rate[1] - rate_exact1)/rate_exact1 ) > tol )
    {
      std::cout << "Error: Mismatch in rate values." << std::endl
		<< "rate(T1)   = " << rate[1] << std::endl
		<< "rate_exact = " << rate_exact1 << std::endl
		<< "difference = " << rate[1] - rate_exact1 << std::endl;

      return_flag = 1;
    }

  std::cout << "Kooij rate: " << kooij_rate << std::endl;

  return return_flag;
}


int main()
{
  int returnval = 0;

  returnval = returnval ||
    vectester (std::valarray<float>(2));
  returnval = returnval ||
    vectester (std::valarray<double>(2));
  returnval = returnval ||
    vectester (std::valarray<long double>(2));
#ifdef ANTIOCH_HAVE_EIGEN
  returnval = returnval ||
    vectester (Eigen::Array2f());
  returnval = returnval ||
    vectester (Eigen::Array2d());
  returnval = returnval ||
    vectester (Eigen::Array<long double, 2, 1>());
#endif
#ifdef ANTIOCH_HAVE_METAPHYSICL
  returnval = returnval ||
    vectester (MetaPhysicL::NumberArray<2, float> (0));
  returnval = returnval ||
    vectester (MetaPhysicL::NumberArray<2, double> (0));
  returnval = returnval ||
    vectester (MetaPhysicL::NumberArray<2, long double> (0));
#endif

  return returnval;
}
