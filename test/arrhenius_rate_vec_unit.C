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

#ifdef ANTIOCH_HAVE_VEXCL
#include "vexcl/vexcl.hpp"
#endif

#include "antioch/eigen_utils_decl.h"
#include "antioch/metaphysicl_utils_decl.h"
#include "antioch/valarray_utils_decl.h"
#include "antioch/vexcl_utils_decl.h"

#include "antioch/arrhenius_rate.h"

#include "antioch/eigen_utils.h"
#include "antioch/metaphysicl_utils.h"
#include "antioch/valarray_utils.h"
#include "antioch/vexcl_utils.h"

#ifdef ANTIOCH_HAVE_GRVY
#include "grvy.h"

GRVY::GRVY_Timer_Class gt;
#endif

#include <cmath>
#include <limits>

template <typename PairScalars>
int vectester(const PairScalars& example, const std::string& testname)
{
  using std::abs;
  using std::exp;

  typedef typename Antioch::value_type<PairScalars>::type Scalar;

  const Scalar Cf = 1.4;
  const Scalar Ea = 5.0;

  Antioch::ArrheniusRate<Scalar> arrhenius_rate(Cf,Ea,1.);

  // Construct from example to avoid resizing issues
  PairScalars T = example;
  for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
    {
      T[2*tuple]   = 1500.1;
      T[2*tuple+1] = 1600.1;
    }
  
  const Scalar rate_exact0 = Cf*exp(-Ea/1500.1);
  const Scalar rate_exact1 = Cf*exp(-Ea/1600.1);
  const Scalar derive_exact0 = Ea/(Scalar(1500.1)*Scalar(1500.1)) * Cf * exp(-Ea/Scalar(1500.1));
  const Scalar derive_exact1 = Ea/(Scalar(1600.1)*Scalar(1600.1)) * Cf * exp(-Ea/Scalar(1600.1));

  int return_flag = 0;

#ifdef ANTIOCH_HAVE_GRVY
  gt.BeginTimer(testname);
#endif

  const PairScalars rate = arrhenius_rate(T);
  const PairScalars deriveRate = arrhenius_rate.derivative(T);

#ifdef ANTIOCH_HAVE_GRVY
  gt.EndTimer(testname);
#endif

  const Scalar tol = std::numeric_limits<Scalar>::epsilon()*10;

  for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
    {
      if( abs( (rate[2*tuple] - rate_exact0)/rate_exact0 ) > tol )
        {
          std::cout << "Error: Mismatch in rate values." << std::endl
		    << "rate(T0)   = " << rate[2*tuple] << std::endl
		    << "rate_exact = " << rate_exact0 << std::endl
		    << "difference = " << rate[2*tuple] - rate_exact0 << std::endl;

          return_flag = 1;
	  break;
        }

      if( abs( (rate[2*tuple+1] - rate_exact1)/rate_exact1 ) > tol )
        {
          std::cout << "Error: Mismatch in rate values." << std::endl
		    << "rate(T1)   = " << rate[2*tuple+1] << std::endl
		    << "rate_exact = " << rate_exact1 << std::endl
		    << "difference = " << rate[2*tuple+1] - rate_exact1 << std::endl;

          return_flag = 1;
	  break;
        }
    }
  if( abs( (deriveRate[0] - derive_exact0)/derive_exact0 ) > tol )
    {
      std::cout << std::scientific << std::setprecision(16)
                << "Error: Mismatch in rate derivative values." << std::endl
		<< "drate_dT(T0) = " << deriveRate[0] << std::endl
		<< "derive_exact = " << derive_exact0 << std::endl;

      return_flag = 1;
    }
  if( abs( (deriveRate[1] - derive_exact1)/derive_exact1 ) > tol )
    {
      std::cout << std::scientific << std::setprecision(16)
                << "Error: Mismatch in rate derivative values." << std::endl
		<< "drate_dT(T1) = " << deriveRate[1] << std::endl
		<< "derive_exact = " << derive_exact1 << std::endl;

      return_flag = 1;
    }

  std::cout << "Arrhenius rate: " << arrhenius_rate << std::endl;

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
  vex::Context ctx (vex::Filter::DoublePrecision);

  returnval = returnval ||
    vectester (vex::vector<float> (ctx, 2*ANTIOCH_N_TUPLES), "vex::vector<float>");
  returnval = returnval ||
    vectester (vex::vector<double> (ctx, 2*ANTIOCH_N_TUPLES), "vex::vector<double>");
#endif

#ifdef ANTIOCH_HAVE_GRVY
  gt.Finalize();
  gt.Summarize();
#endif

  return returnval;
}
