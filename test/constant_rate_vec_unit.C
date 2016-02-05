
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

#include "antioch/constant_rate.h"

#include "antioch/eigen_utils.h"
#include "antioch/metaphysicl_utils.h"
#include "antioch/valarray_utils.h"

#include <cmath>
#include <limits>

#ifdef ANTIOCH_HAVE_GRVY
#include "grvy.h"

GRVY::GRVY_Timer_Class gt;
#endif

template <typename PairScalars>
int check_rate_and_derivative(const PairScalars & rate_exact, const PairScalars & derive_exact,
                              const PairScalars & rate,       const PairScalars & derive, const PairScalars & T)
{
    typedef typename Antioch::value_type<PairScalars>::type  Scalar;
    const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 2;

    int return_flag(0);
   for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
   {
    if( abs( (rate[2*tuple] - rate_exact[2*tuple])/rate_exact[2*tuple] ) > tol )
      {
        std::cout << std::scientific << std::setprecision(16)
                  << "Error: Mismatch in rate values." << std::endl
                  << "T = " << T[2*tuple] << " K" << std::endl
                  << "rate(T) = " << rate[2*tuple] << std::endl
                  << "rate_exact = " << rate_exact[2*tuple] << std::endl
                  << "relative difference = " <<  abs( (rate[2*tuple] - rate_exact[2*tuple])/rate_exact[2*tuple] ) << std::endl
                  << "tolerance = " <<  tol << std::endl;

        return_flag = 1;
      }
    if( abs( (rate[2*tuple+1] - rate_exact[2*tuple+1])/rate_exact[2*tuple+1] ) > tol )
      {
        std::cout << std::scientific << std::setprecision(16)
                  << "Error: Mismatch in rate values." << std::endl
                  << "T = " << T[2*tuple+1] << " K" << std::endl
                  << "rate(T) = " << rate[2*tuple] << std::endl
                  << "rate_exact = " << rate_exact[2*tuple+1] << std::endl
                  << "relative difference = " <<  abs( (rate[2*tuple] - rate_exact[2*tuple+1])/rate_exact[2*tuple+1] ) << std::endl
                  << "tolerance = " <<  tol << std::endl;

        return_flag = 1;
      }
    if( abs( (derive[2*tuple] - derive_exact[2*tuple])/derive_exact[2*tuple] ) > tol )
      {
        std::cout << std::scientific << std::setprecision(16)
                  << "Error: Mismatch in rate derivative values." << std::endl
                  << "T = " << T[2*tuple] << " K" << std::endl
                  << "drate_dT(T) = " << derive[2*tuple] << std::endl
                  << "derive_exact = " << derive_exact[2*tuple] << std::endl
                  << "relative difference = " <<  abs( (derive[2*tuple] - derive_exact[2*tuple])/derive_exact[2*tuple] ) << std::endl
                  << "tolerance = " <<  tol << std::endl;

        return_flag = 1;
      }
    if( abs( (derive[2*tuple+1] - derive_exact[2*tuple+1])/derive_exact[2*tuple+1] ) > tol )
      {
        std::cout << std::scientific << std::setprecision(16)
                  << "Error: Mismatch in rate derivative values." << std::endl
                  << "T = " << T[2*tuple+1] << " K" << std::endl
                  << "drate_dT(T) = " << derive[2*tuple+1] << std::endl
                  << "derive_exact = " << derive_exact[2*tuple+1] << std::endl
                  << "relative difference = " <<  abs( (derive[2*tuple+1] - derive_exact[2*tuple+1])/derive_exact[2*tuple+1] ) << std::endl
                  << "tolerance = " <<  tol << std::endl;

        return_flag = 1;
      }

   }
   return return_flag;
}

template <typename PairScalars>
int vectester(const PairScalars& example, const std::string & testname)
{
  using std::abs;

  typedef typename Antioch::value_type<PairScalars>::type Scalar;

  const Scalar Cf = 1.4L;

  Antioch::ConstantRate<Scalar> constant_rate(Cf);

  // Construct from example to avoid resizing issues
  PairScalars T = example;
  PairScalars rate_exact = example;
  PairScalars derive_exact =  Antioch::zero_clone(example);
  for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
  {
    T[2*tuple] = 1500.1L;
    T[2*tuple+1] = 1600.1L;
    rate_exact[2*tuple] = Cf;
    rate_exact[2*tuple+1] = Cf;
  }
  Antioch::KineticsConditions<PairScalars> cond(T);

  int return_flag = 0;

#ifdef ANTIOCH_HAVE_GRVY
  gt.BeginTimer(testname);
#endif

// KineticsConditions
  PairScalars rate = constant_rate(cond);//Antioch::zero_clone(T.T());
  PairScalars derive = constant_rate.derivative(cond);//Antioch::zero_clone(T.T());

#ifdef ANTIOCH_HAVE_GRVY
  gt.EndTimer(testname);
#endif

  return_flag = check_rate_and_derivative(rate_exact, derive_exact, rate, derive,T) || return_flag;

#ifdef ANTIOCH_HAVE_GRVY
  gt.BeginTimer(testname);
#endif

  constant_rate.rate_and_derivative(cond,rate,derive);

#ifdef ANTIOCH_HAVE_GRVY
  gt.EndTimer(testname);
#endif

  return_flag = check_rate_and_derivative(rate_exact, derive_exact, rate, derive,T) || return_flag;

// T

#ifdef ANTIOCH_HAVE_GRVY
  gt.BeginTimer(testname);
#endif

  rate = constant_rate(T);
  derive = constant_rate.derivative(T);

#ifdef ANTIOCH_HAVE_GRVY
  gt.EndTimer(testname);
#endif

  return_flag = check_rate_and_derivative(rate_exact, derive_exact, rate, derive,T) || return_flag;

#ifdef ANTIOCH_HAVE_GRVY
  gt.BeginTimer(testname);
#endif

  constant_rate.rate_and_derivative(T,rate,derive);

#ifdef ANTIOCH_HAVE_GRVY
  gt.EndTimer(testname);
#endif

  return_flag = check_rate_and_derivative(rate_exact, derive_exact, rate, derive,T) || return_flag;

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
