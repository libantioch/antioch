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
#include "antioch/vector_utils_decl.h"

#include "antioch/pure_species_viscosity.h"
#include "antioch/gsl_spliner.h"

#include "antioch/eigen_utils.h"
#include "antioch/metaphysicl_utils.h"
#include "antioch/valarray_utils.h"
#include "antioch/vexcl_utils.h"
#include "antioch/vector_utils.h"

#include "antioch/physical_constants.h"

#ifdef ANTIOCH_HAVE_GRVY
#include "grvy.h"

GRVY::GRVY_Timer_Class gt;
#endif

#include <cmath>
#include <limits>

template <typename Scalar, typename Element>
int test_viscosity( const Element & mu, const Scalar  & mu_exact, const Scalar & tol )
{
  using std::abs;

  int return_flag = 0;

  const double rel_error = abs( (mu - mu_exact)/mu_exact);

  if( rel_error  > tol )
    {
      std::cerr << std::setprecision(15) << std::scientific;
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
  using std::abs;
  using std::exp;

  typedef typename Antioch::value_type<PairScalars>::type Scalar;

// value for N2
  const Scalar LJ_depth(97.530L);
  const Scalar LJ_diameter(3.621L);
  const Scalar dipole_moment(0.L);
  const Scalar mass(28.016e-3L/Antioch::Constants::Avogadro<Scalar>());

  Antioch::PureSpeciesViscosity<Scalar,Antioch::GSLSpliner> mu_sp(LJ_depth,LJ_diameter,dipole_moment,mass);

  // Construct from example to avoid resizing issues
  PairScalars T = example;
  for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
    {
      T[2*tuple]   = 1500.1;
      T[2*tuple+1] = 1600.1;
    }
  
  const Scalar mu_exact0 = 0.0000417395098853601937871105407365424874568203066945573066;
  const Scalar mu_exact1 = 0.0000431082906407300464864929380770860406892212065614947462;

  int return_flag = 0;

#ifdef ANTIOCH_HAVE_GRVY
  gt.BeginTimer(testname);
#endif

  PairScalars mu = mu_sp(T) * mu_sp.Stockmayer(T);

#ifdef ANTIOCH_HAVE_GRVY
  gt.EndTimer(testname);
#endif

  const Scalar tol = (std::numeric_limits<Scalar>::epsilon()*10 < 5e-17)?5e-17:
                       std::numeric_limits<Scalar>::epsilon()*10;

  for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
    {
          return_flag = test_viscosity( mu[2*tuple], mu_exact0, tol ) || return_flag;
          return_flag = test_viscosity( mu[2*tuple+1], mu_exact1, tol ) || return_flag;
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
