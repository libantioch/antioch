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
#include "antioch/vector_utils_decl.h"

#include "antioch/molecular_binary_diffusion.h"
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

template <typename Scalar, typename Element, typename ElementOrScalar>
int test_diff( const Element & dij, const ElementOrScalar & dij_exact, const Scalar & tol, const std::string & words )
{
  using std::abs;

  int return_flag = 0;

  const double rel_error = abs( (dij - dij_exact)/dij_exact);

  if( rel_error  > tol )
    {
      std::cerr << std::setprecision(15) << std::scientific;
      std::cerr << "Error: Mismatch in bimolecular coefficient of " << words << std::endl
		<< "Dij(T)    = " << dij << std::endl
		<< "Dij_exact = " << dij_exact << std::endl
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

//  N2
   const Scalar N2_LJ_eps(97.530L);
   const Scalar N2_LJ_depth(3.621L);
   const Scalar N2_dipole(0.L);
   const Scalar N2_polar(1.760L);
   const Scalar N2_Zrot(4.L);
   const Scalar N2_mass(28.016e-3L);
// CH4
   const Scalar CH4_LJ_eps(141.400L);
   const Scalar CH4_LJ_depth(3.746L);
   const Scalar CH4_dipole(0.L);
   const Scalar CH4_polar(2.6L);
   const Scalar CH4_Zrot(13.L);
   const Scalar CH4_mass(16.043e-3L);
// H2O
   const Scalar H2O_LJ_eps(572.4L);
   const Scalar H2O_LJ_depth(2.605L);
   const Scalar H2O_dipole(1.844L);
   const Scalar H2O_polar(0.L);
   const Scalar H2O_Zrot(4.L);
   const Scalar H2O_mass(18.016e-3L);

  Antioch::TransportSpecies<Scalar> N2(0,N2_LJ_eps,N2_LJ_depth,N2_dipole,N2_polar,N2_Zrot,N2_mass),
                                    CH4(1,CH4_LJ_eps,CH4_LJ_depth,CH4_dipole,CH4_polar,CH4_Zrot,CH4_mass),
                                    H2O(2,H2O_LJ_eps,H2O_LJ_depth,H2O_dipole,H2O_polar,H2O_Zrot,H2O_mass);


   Antioch::MolecularBinaryDiffusion<Scalar,Antioch::GSLSpliner> 
                D00(N2,N2),  D01(N2,CH4),  D02(N2,H2O),
                D10(CH4,N2), D11(CH4,CH4), D12(CH4,H2O),
                D20(H2O,N2), D21(H2O,CH4), D22(H2O,H2O);


  // Construct from example to avoid resizing issues
  PairScalars T = example;
  PairScalars cTot = example;
  for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
    {
      T[2*tuple]      = 1500.1;
      T[2*tuple+1]    = 1600.1;
      cTot[2*tuple]   = 5e-7;
      cTot[2*tuple+1] = 5e-7;
    }
  // bc gives
  const Scalar D00_exact0 = 3.575629059282712203350417538824313603e3;
  const Scalar D11_exact0 = 4.415035849326582722446637772820322872e3;
  const Scalar D22_exact0 = 8.615250909281137767894964155009068175e3;
  const Scalar D10_exact0 = 4.048999169845892614961423528844221310e3;
  const Scalar D20_exact0 = 5.468107000169991297723144486054211050e3;
  const Scalar D12_exact0 = 5.973341001459783642751059656941311432e3;

  const Scalar D00_exact1 = 3.692886119994007408894312228191265622e3;
  const Scalar D11_exact1 = 4.559819918938905357528221639097498410e3;
  const Scalar D22_exact1 = 8.897774342826358536851687124754748934e3;
  const Scalar D10_exact1 = 4.181779649478152213931731689698388959e3;
  const Scalar D20_exact1 = 5.647424861129375458731724279131726416e3;
  const Scalar D12_exact1 = 6.169227206892386749039436637869885641e3;

#ifdef ANTIOCH_HAVE_GRVY
  gt.BeginTimer(testname);
#endif
  
  PairScalars d00 = D00(T,cTot) * D00.Stockmayer(T);
  PairScalars d11 = D11(T,cTot) * D11.Stockmayer(T);
  PairScalars d22 = D22(T,cTot) * D22.Stockmayer(T);

  PairScalars d10 = D10(T,cTot) * D10.Stockmayer(T);
  PairScalars d01 = D01(T,cTot) * D01.Stockmayer(T);
  PairScalars d20 = D20(T,cTot) * D20.Stockmayer(T);
  PairScalars d02 = D02(T,cTot) * D02.Stockmayer(T);
  PairScalars d21 = D21(T,cTot) * D21.Stockmayer(T);
  PairScalars d12 = D12(T,cTot) * D12.Stockmayer(T);

  int return_flag = 0;

#ifdef ANTIOCH_HAVE_GRVY
  gt.EndTimer(testname);
#endif

  const Scalar tol = (std::numeric_limits<Scalar>::epsilon()*10 < 5e-16)?5e-16:
                       std::numeric_limits<Scalar>::epsilon()*10;

  for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
    {
        return_flag = test_diff( d10[2*tuple] , d01[2*tuple], tol, "N2 - CH4 symetry" )   || return_flag;
        return_flag = test_diff( d20[2*tuple] , d02[2*tuple], tol, "N2 - H2O symetry" )   || return_flag;
        return_flag = test_diff( d12[2*tuple] , d21[2*tuple], tol, "H2O - CH4 symetry" )   || return_flag;
//
        return_flag = test_diff( d00[2*tuple] , D00_exact0, tol, "N2 - N2" )   || return_flag;
        return_flag = test_diff( d11[2*tuple] , D11_exact0, tol, "CH4 - CH4" ) || return_flag;
        return_flag = test_diff( d22[2*tuple] , D22_exact0, tol, "H2O - H2O" ) || return_flag;
        return_flag = test_diff( d10[2*tuple] , D10_exact0, tol, "N2 - CH4" )  || return_flag;
        return_flag = test_diff( d12[2*tuple] , D12_exact0, tol, "CH4 - H2O" ) || return_flag;
        return_flag = test_diff( d20[2*tuple] , D20_exact0, tol, "N2 - H2O" )  || return_flag;
//
        return_flag = test_diff( d00[2*tuple+1] , D00_exact1, tol, "N2 - N2" )   || return_flag;
        return_flag = test_diff( d11[2*tuple+1] , D11_exact1, tol, "CH4 - CH4" ) || return_flag;
        return_flag = test_diff( d22[2*tuple+1] , D22_exact1, tol, "H2O - H2O" ) || return_flag;
        return_flag = test_diff( d10[2*tuple+1] , D10_exact1, tol, "N2 - CH4" )  || return_flag;
        return_flag = test_diff( d12[2*tuple+1] , D12_exact1, tol, "CH4 - H2O" ) || return_flag;
        return_flag = test_diff( d20[2*tuple+1] , D20_exact1, tol, "N2 - H2O" )  || return_flag;
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
