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

#include "antioch/stat_mech_thermo.h"
#include "antioch/chemical_mixture.h"
#include "antioch/pure_species_thermal_conductivity.h"

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

template <typename Scalar, typename Element>
int test_k( const Element & k, const Scalar & k_exact, const Scalar & tol )
{
  using std::abs;

  int return_flag = 0;

  const Scalar rel_error = abs( (k - k_exact)/k_exact);

  if( rel_error  > tol )
    {
      std::cerr << std::setprecision(20) << std::scientific
                << "Error: Mismatch in thermal conductivity" << std::endl
		<< "k       = " << k << std::endl
		<< "k_exact = " << k_exact << std::endl
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

  std::vector<std::string> species_str_list;
  const unsigned int n_species = 1;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );

  const Scalar LJ_depth_N2 = 97.53L;
  const Scalar Z_298 = 4.0L;

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );

  Antioch::StatMechThermodynamics<Scalar> thermo( chem_mixture );

  Antioch::PureSpeciesThermalConductivity<Antioch::StatMechThermodynamics<Scalar>, Scalar > k( thermo, Z_298, LJ_depth_N2);

  // Construct from example to avoid resizing issues
  PairScalars mu  = example;
  PairScalars dss = example;
  PairScalars rho = example;
  PairScalars T = example;
  for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
    {
      T[2*tuple]     = 1500.1;
      T[2*tuple+1]   = 1600.1;
      mu[2*tuple]    = 3.14e-3;
      mu[2*tuple+1]  = 3.14e-3;
      dss[2*tuple]   = 5.23e-5;
      dss[2*tuple+1] = 5.23e-5;
      rho[2*tuple]   = 1.4;
      rho[2*tuple+1] = 1.4;
    }
  
  const Scalar k_exact0 = 3.19434291925996020233442116364270671873509961381739235964650684;
  const Scalar k_exact1 = 3.2021908709042895344166123247634817043607657971633006497668;

  int return_flag = 0;

#ifdef ANTIOCH_HAVE_GRVY
  gt.BeginTimer(testname);
#endif

  const PairScalars k_ps = k(0,mu,T,rho,dss);

#ifdef ANTIOCH_HAVE_GRVY
  gt.EndTimer(testname);
#endif

  const Scalar tol = (std::numeric_limits<Scalar>::epsilon()*10 < 7e-17)?
                        7e-17:
                        std::numeric_limits<Scalar>::epsilon()*10;

  for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
    {

      return_flag = test_k( k_ps[2*tuple], k_exact0, tol ) || return_flag;
      return_flag = test_k( k_ps[2*tuple+1], k_exact1, tol ) || return_flag;
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
