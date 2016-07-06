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

#include "antioch/vanthoff_rate.h"

#include "antioch/eigen_utils.h"
#include "antioch/metaphysicl_utils.h"
#include "antioch/valarray_utils.h"
#include "antioch/vexcl_utils.h"

#include <cmath>
#include <limits>

#ifdef ANTIOCH_HAVE_GRVY
#include "grvy.h"

GRVY::GRVY_Timer_Class gt;
#endif

template <typename Scalar>
bool excessive_difference(Scalar T, Scalar approx, Scalar exact, Scalar tol, const std::string & name)
{
  using std::abs;
  using std::max;

  Scalar relative_difference = abs(approx - exact)/max(abs(exact),tol);
  if (relative_difference > tol)
    std::cout << std::scientific << std::setprecision(16)
             << "Error: Mismatch in " << name << '.' << std::endl
             << "T = " << T << " K" << std::endl
             << "approx(T) = " << approx << std::endl
             << "exact(T)  = " << exact << std::endl
             << "relative difference = " << relative_difference << std::endl
             << "tolerance = " <<  tol << std::endl;

  return (relative_difference > tol);
}



template <typename PairScalars>
int check_rate_and_derivative(const PairScalars & rate_exact, const PairScalars & derive_exact,
                              const PairScalars & rate,       const PairScalars & derive, const PairScalars & T)
{
    typedef typename Antioch::value_type<PairScalars>::type  Scalar;
    const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 2;

    int return_flag(0);
   for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
   {
     if( excessive_difference<Scalar>(T[2*tuple],   rate[2*tuple],     rate_exact[2*tuple],     tol, "rate values") )
       return_flag = 1;
     if( excessive_difference<Scalar>(T[2*tuple+1], rate[2*tuple+1],   rate_exact[2*tuple+1],   tol, "rate values") )
       return_flag = 1;
     if( excessive_difference<Scalar>(T[2*tuple],   derive[2*tuple],   derive_exact[2*tuple],   tol, "derivative values") )
       return_flag = 1;
     if( excessive_difference<Scalar>(T[2*tuple+1], derive[2*tuple+1], derive_exact[2*tuple+1], tol, "derivative values") )
       return_flag = 1;
   }
   return return_flag;
}

template <typename PairScalars>
int vectester(const PairScalars& example, const std::string & testname)
{
  using std::abs;
  using std::exp;
  using std::pow;

  typedef typename Antioch::value_type<PairScalars>::type Scalar;

  const Scalar Cf = 1.4;
  const Scalar eta = 1.2;
  const Scalar Ea = 5.0;
  const Scalar D = -2.5;

  Antioch::VantHoffRate<Scalar> vanthoff_rate(Cf,eta,Ea,D,1.,1.);

  // Construct from example to avoid resizing issues
  PairScalars T = example;
  PairScalars rate_exact = example;
  PairScalars derive_exact = example;
  for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
    {
      T[2*tuple] = 1500.1;
      T[2*tuple+1] = 1600.1;
      rate_exact[2*tuple] = Cf*pow(Scalar(1500.1),eta)*exp(-Ea/1500.1+D*1500.1);
      rate_exact[2*tuple+1] = Cf*pow(Scalar(1600.1),eta)*exp(-Ea/1600.1+D*1600.1);
      derive_exact[2*tuple] = Cf * pow(Scalar(1500.1),eta) * exp(-Ea/Scalar(1500.1) + D*Scalar(1500.1)) * (D + eta/Scalar(1500.1) + Ea/(Scalar(1500.1)*Scalar(1500.1)));
      derive_exact[2*tuple+1] = Cf * pow(Scalar(1600.1),eta) * exp(-Ea/Scalar(1600.1) + D*Scalar(1600.1)) * (D + eta/Scalar(1600.1) + Ea/(Scalar(1600.1)*Scalar(1600.1)));
    }
  Antioch::KineticsConditions<PairScalars> cond(T);

  int return_flag = 0;

// KineticsConditions
#ifdef ANTIOCH_HAVE_GRVY
  gt.BeginTimer(testname);
#endif

  PairScalars rate = vanthoff_rate(cond);
  PairScalars derive = vanthoff_rate.derivative(cond);

#ifdef ANTIOCH_HAVE_GRVY
  gt.EndTimer(testname);
#endif

  return_flag = check_rate_and_derivative(rate_exact, derive_exact, rate, derive,T) || return_flag;

#ifdef ANTIOCH_HAVE_GRVY
  gt.BeginTimer(testname);
#endif

  vanthoff_rate.rate_and_derivative(cond,rate,derive);

#ifdef ANTIOCH_HAVE_GRVY
  gt.EndTimer(testname);
#endif

  return_flag = check_rate_and_derivative(rate_exact, derive_exact, rate, derive,T) || return_flag;

// T
#ifdef ANTIOCH_HAVE_GRVY
  gt.BeginTimer(testname);
#endif

  rate = vanthoff_rate(T);
  derive = vanthoff_rate.derivative(T);

#ifdef ANTIOCH_HAVE_GRVY
  gt.EndTimer(testname);
#endif

  return_flag = check_rate_and_derivative(rate_exact, derive_exact, rate, derive,T) || return_flag;

#ifdef ANTIOCH_HAVE_GRVY
  gt.BeginTimer(testname);
#endif

  vanthoff_rate.rate_and_derivative(cond,rate,derive);

#ifdef ANTIOCH_HAVE_GRVY
  gt.EndTimer(testname);
#endif

  return_flag = check_rate_and_derivative(rate_exact, derive_exact, rate, derive,T) || return_flag;

  std::cout << "Van't Hoff rate: " << vanthoff_rate << std::endl;

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
