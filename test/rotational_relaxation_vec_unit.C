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
#include "antioch/vexcl_utils_decl.h"

#include "antioch/rotational_relaxation.h"

#include "antioch/eigen_utils.h"
#include "antioch/metaphysicl_utils.h"
#include "antioch/valarray_utils.h"
#include "antioch/vexcl_utils.h"

#include "antioch/cmath_shims.h"

#ifdef ANTIOCH_HAVE_GRVY
#include "grvy.h"

GRVY::GRVY_Timer_Class gt;
#endif

#include <cmath>
#include <limits>

const long double pi(3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825);

template <typename Scalar>
Scalar F(const Scalar & x)
{
  return 1.L + Antioch::ant_pow(pi,(float)1.5)/2.L  * Antioch::ant_sqrt(x)
             + (pi*pi/4.L + 2.L) * x
             + Antioch::ant_pow(pi * x,(float)1.5);
}

template <typename Scalar>
Scalar z(const Scalar & T, const Scalar & eps_kb, const Scalar & z_298)
{
std::cout << "iiiii" << std::endl;
  return z_298 * F(eps_kb / 298.L) / F(eps_kb / T);
}


template <typename PairScalars>
int vectester(const PairScalars& example, const std::string& testname)
{

std::cout << "entering " << testname << std::endl;
  typedef typename Antioch::value_type<PairScalars>::type Scalar;

  const Scalar eps_kb = 97.53L; // N2 value
  const Scalar z_298 = 4.0L; // N2 value

std::cout << "0" << std::endl;
  Antioch::RotationalRelaxation<Scalar> rot(z_298,eps_kb);

std::cout << "1" << std::endl;
  // Construct from example to avoid resizing issues
  PairScalars T = example;
  for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
    {
      T[2*tuple]   = 1500.1L;
      T[2*tuple+1] = 1600.1L;
    }
  
  int return_flag = 0;

std::cout << "2" << std::endl;
#ifdef ANTIOCH_HAVE_GRVY
  gt.BeginTimer(testname);
#endif

std::cout << "3" << std::endl;
  const PairScalars Z = rot(T);
std::cout << "4" << std::endl;

#ifdef ANTIOCH_HAVE_GRVY
  gt.EndTimer(testname);
#endif

  const Scalar tol = (std::numeric_limits<Scalar>::epsilon()*10 < 5e-17)?
                        5e-17:
                        std::numeric_limits<Scalar>::epsilon()*10;

  for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
    {
std::cout << "mmmmm " << tuple << " / " << ANTIOCH_N_TUPLES << "  " << testname << std::endl;
      const Scalar Z_exact0 = z( (Scalar)1500.1, eps_kb, z_298);
std::cout << "ppppp" << std::endl;
      const Scalar Z_exact1 = z( (Scalar)1600.1, eps_kb, z_298);

      if( abs( (Z[2*tuple] - Z_exact0)/Z_exact0 ) > tol )
        {
          std::cout << std::scientific << std::setprecision(20)
                    << "Error: Mismatch in Z values in test " << testname << std::endl
		    << "Z(T0)   = "             << Z[2*tuple]                                 << std::endl
		    << "Z_exact = "             << Z_exact0                                   << std::endl
		    << "absolute difference = " << Z[2*tuple] - Z_exact0                      << std::endl
		    << "relative difference = " << std::abs(Z[2*tuple] - Z_exact0) / Z_exact0 << std::endl
		    << "tolerance = "           << tol << std::endl;

          return_flag = 1;
	  break;
        }

      if( abs( (Z[2*tuple+1] - Z_exact1)/Z_exact1 ) > tol )
        {
          std::cout << std::scientific << std::setprecision(20)
                    << "Error: Mismatch in Z values in test " << testname << std::endl
		    << "Z(T1)   = "             << Z[2*tuple+1]                                 << std::endl
		    << "Z_exact = "             << Z_exact1                                     << std::endl
		    << "absolute difference = " << Z[2*tuple+1] - Z_exact1                      << std::endl
		    << "relative difference = " << std::abs(Z[2*tuple+1] - Z_exact1) / Z_exact1 << std::endl
		    << "tolerance = "           << tol << std::endl;

          return_flag = 1;
	  break;
        }
    }

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
