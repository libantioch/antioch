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
#include "antioch/vector_utils_decl.h"

#include "antioch/gsl_spliner.h"

#include "antioch/eigen_utils.h"
#include "antioch/metaphysicl_utils.h"
#include "antioch/valarray_utils.h"
#include "antioch/vexcl_utils.h"
#include "antioch/vector_utils.h"

#ifdef ANTIOCH_HAVE_GRVY
#include "grvy.h"

GRVY::GRVY_Timer_Class gt;
#endif

#include <cmath>
#include <limits>

template <typename Scalar,typename Element>
int check_value(const Element & ref, const Element & candidate, const Element & x, const std::string & words)
{
 /*because not the real test function impossible due to boundary conditions
  */
  const Scalar tol = 5e-3;//std::numeric_limits<Scalar>::epsilon() * 10.;

  if(std::abs((ref - candidate)/ref) > tol)
  {
      std::cerr << std::scientific << std::setprecision(15);
      std::cerr << "ERROR in gsl spline test at point: " << x << "; " << words << std::endl
                << "  reference = " << ref << std::endl
                << "  candidate = " << candidate << std::endl
                << "  relative difference = " << std::abs((ref - candidate) / ref) << std::endl
                << "  absolute difference = " << std::abs(ref - candidate) << std::endl
                << "  tolerance = " << tol << std::endl;
      return 1;
  }

  return 0;
}

template <typename Scalar>
Scalar function(const Scalar x)
{
  Scalar ten  = Antioch::constant_clone(x,10);
  Scalar two  = Antioch::constant_clone(x,2);
  Scalar five = Antioch::constant_clone(x,5);

  return ten + five * x + ten * x * x - two * x * x * x;
}

template <typename Scalar>
void fill_ref(std::vector<Scalar> & x_ref, std::vector<Scalar> & y_ref,
              unsigned int n_data, const Scalar & min,
              const Scalar & max)
{
  for(unsigned int i = 0; i < n_data; i++)
  {
     x_ref[i] = min + (Scalar)(i) * (max - min) / (Scalar)(n_data-1);
     y_ref[i] = function(x_ref[i]);
  }
}

template <typename PairScalars>
int vectester(const PairScalars& example, const std::string& testname)
{
  typedef typename Antioch::value_type<PairScalars>::type Scalar;

  const unsigned int n_data(40);

  std::vector<Scalar> x_ref(n_data,0),y_ref(n_data,0);

  const Scalar min = -5L;
  const Scalar max = 8L;
  fill_ref(x_ref,y_ref,n_data,min, max);

  Antioch::GSLSpliner gsl_spline(x_ref,y_ref);

  // Construct from example to avoid resizing issues
  PairScalars x = example;
  for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
    {
      x[2*tuple]   = -3.5;
      x[2*tuple+1] = 5.1;
    }
  
  int return_flag = 0;

#ifdef ANTIOCH_HAVE_GRVY
  gt.BeginTimer(testname);
#endif

  const PairScalars gsl = gsl_spline.interpolated_value(x);

#ifdef ANTIOCH_HAVE_GRVY
  gt.EndTimer(testname);
#endif

    const PairScalars exact = function(x);
    for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
      {
          return_flag = check_value<Scalar>(exact[2*tuple],   gsl[2*tuple],   x[2*tuple],   "gsl vectorized") || return_flag;
          return_flag = check_value<Scalar>(exact[2*tuple+1], gsl[2*tuple+1], x[2*tuple+1], "gsl vectorized") || return_flag;
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
