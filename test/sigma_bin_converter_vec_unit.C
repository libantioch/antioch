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

#include "antioch/vector_utils_decl.h"
#include "antioch/eigen_utils_decl.h"
#include "antioch/metaphysicl_utils_decl.h"
#include "antioch/valarray_utils_decl.h"
#include "antioch/vexcl_utils_decl.h"

#include "antioch/sigma_bin_converter.h"
#include "antioch/antioch_asserts.h"
#include "antioch/physical_constants.h"

#include "antioch/vector_utils.h"
#include "antioch/eigen_utils.h"
#include "antioch/metaphysicl_utils.h"
#include "antioch/valarray_utils.h"
#include "antioch/vexcl_utils.h"

#ifdef ANTIOCH_HAVE_GRVY
#include "grvy.h"

GRVY::GRVY_Timer_Class gt;
#endif

// C++
#include <cmath>
#include <limits>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

/*
  x[0] = 1.;  y[0] = 1.;
  x[1] = 2.;  y[1] = 2.;
  x[2] = 3.;  y[2] = 5.;
  x[3] = 4.;  y[3] = 8.;
  x[4] = 5.;  y[4] = 6.;
  x[5] = 6.;  y[5] = 10.;
  x[6] = 7.;  y[6] = 7.;
  x[7] = 8.;  y[7] = 4.;
  x[8] = 9.;  y[8] = 2.;
  x[9] = 10.; y[9] = 0.3; // right stairs, this value is useless
*/
template <typename PairScalars, typename VectorPairScalar>
void make_custom(unsigned int choice, const PairScalars & ex, VectorPairScalar & x, VectorPairScalar & y)
{
  // sum_{bin} Delta x * y
  switch(choice)
  {
      case(0): // 9 bin contained in ref -> [2.5;8.5] within [1;10]
      {
       x.resize(9,ex);   y.resize(9,ex);
       for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
       {
          x[0][2*tuple] = 2.50L;  y[0][2*tuple] = 2.25L/0.75L; // 0.50][2  + 0.25][5
          x[1][2*tuple] = 3.25L;  y[1][2*tuple] = 3.75L/0.75L; // 0.75][5
          x[2][2*tuple] = 4.00L;  y[2][2*tuple] = 6.00L/0.75L; // 0.75][8
          x[3][2*tuple] = 4.75L;  y[3][2*tuple] = 5.00L/0.75L; // 0.25][8  + 0.5 ][6
          x[4][2*tuple] = 5.50L;  y[4][2*tuple] = 5.50L/0.75L; // 0.25][6  + 0.25][10
          x[5][2*tuple] = 6.25L;  y[5][2*tuple] = 7.50L/0.75L; // 0.75][10
          x[6][2*tuple] = 7.00L;  y[6][2*tuple] = 5.25L/0.75L; // 0.75][7
          x[7][2*tuple] = 7.75L;  y[7][2*tuple] = 3.75L/0.75L; // 0.25][7 + 0.5][4
          x[8][2*tuple] = 8.50L;  y[8][2*tuple] = 0.00L/0.75L; // 0. (right stairs)

      // 9 bin outside ref -> [0;12] containing [1;10]
          
          x[0][2*tuple + 1] = 0.00L;  y[0][2*tuple + 1] = 0.50L/1.50L; // 0.5][1
          x[1][2*tuple + 1] = 1.50L;  y[1][2*tuple + 1] = 2.50L/1.50L; // 0.5][1  + 2.0][1
          x[2][2*tuple + 1] = 3.00L;  y[2][2*tuple + 1] = 9.00L/1.50L; // 1.0][5  + 0.5][8
          x[3][2*tuple + 1] = 4.50L;  y[3][2*tuple + 1] = 10.0L/1.50L; // 0.5][8  + 1.0][6
          x[4][2*tuple + 1] = 6.00L;  y[4][2*tuple + 1] = 13.5L/1.50L; // 1.0][10 + 0.5][7
          x[5][2*tuple + 1] = 7.50L;  y[5][2*tuple + 1] = 7.50L/1.50L; // 0.50][7 + 1.0][4
          x[6][2*tuple + 1] = 9.00L;  y[6][2*tuple + 1] = 2.00L/1.50L; // 1.0][2 
          x[7][2*tuple + 1] = 10.5L;  y[7][2*tuple + 1] = 0.00L/1.50L; // 0
          x[8][2*tuple + 1] = 12.0L;  y[8][2*tuple + 1] = 0.00L/1.50L; // 0. (right stairs)
       }
       break;
      }
      case(1): // 5 bins beyond only min -> [-1;3.8] in [0;10]
      {
       x.resize(5,ex);    y.resize(5,ex);
       for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
       {
          x[0][2*tuple] = -1.00L;  y[0][2*tuple] = 0.00L/1.20L; // 0
          x[1][2*tuple] =  0.20L;  y[1][2*tuple] = 0.40L/1.20L; // 0.4][1
          x[2][2*tuple] =  1.40L;  y[2][2*tuple] = 1.80L/1.20L; // 0.6][1  + 0.6][2
          x[3][2*tuple] =  2.60L;  y[3][2*tuple] = 4.80L/1.20L; // 0.4][2  + 0.8][5
          x[4][2*tuple] =  3.80L;  y[4][2*tuple] = 00.0L/1.20L; // 0. (right stairs)
       
        // 6 bins beyond only max -> [2;10.75] in [0;10]
          x[0][2*tuple + 1] =  2.0000L;  y[0][2*tuple + 1] =  8.50L/2.1875L; // 1.0000][2  + 1.0000][5 + 0.1875][8
          x[1][2*tuple + 1] =  4.1875L;  y[1][2*tuple + 1] = 16.25L/2.1875L; // 0.8125][8  + 1.0000][6 + 0.3750][10
          x[2][2*tuple + 1] =  6.3750L;  y[2][2*tuple + 1] = 15.50L/2.1875L; // 0.6250][10 + 1.0000][7 + 0.5625][4
          x[3][2*tuple + 1] =  8.5625L;  y[3][2*tuple + 1] =  3.75L/2.1875L; // 0.4375][4  + 1.0000][2 + 0
          x[4][2*tuple + 1] = 10.7500L;  y[4][2*tuple + 1] =  0.00L/2.1875L; // 0. (right stairs)
        }
        break;
      }
  }

}

template <typename VectorScalar>
void make_reference(VectorScalar & x, VectorScalar & y)
{
  x.resize(10,0);
  y.resize(10,0);
  for(unsigned int i = 0; i < 10; i++)
  {
     x[i] = ((typename Antioch::value_type<VectorScalar>::type)(i+1));
  }
  y[0] = 1.L;
  y[1] = 2.L;
  y[2] = 5.L;
  y[3] = 8.L;
  y[4] = 6.L;
  y[5] = 10.L;
  y[6] = 7.L;
  y[7] = 4.L;
  y[8] = 2.L;
  y[9] = 0.3L;

}


template <typename PairScalars>
int vectester(const PairScalars& example, const std::string& testname)
{
  typedef typename Antioch::value_type<PairScalars>::type Scalar;

  std::vector<Scalar> bin_ref_x,bin_ref_y;

  make_reference(bin_ref_x,bin_ref_y);

  Antioch::SigmaBinConverter<std::vector<Scalar> > bin;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 10;

  int return_flag = 0;

  // 2 * 2  cases:
  //   - custom inside ref
  //   - ref inside custom
  //   - custom beyond only min ref
  //   - custom beyond only max ref
  for(unsigned int i = 0; i < 2; i++)
  {
    std::vector<PairScalars> bin_custom_x, exact_sol_y;
    std::vector<PairScalars> bin_custom_y;
    make_custom(i,example,bin_custom_x,exact_sol_y);

    bin.y_on_custom_grid(bin_ref_x,bin_ref_y,
                         bin_custom_x,bin_custom_y);


    for(unsigned int il = 0; il < bin_custom_x.size() - 1; il++)
    {
      for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
      {
//tuple

        Scalar dist = Antioch::if_else(exact_sol_y[il][2*tuple] < tol,
                                                std::abs(bin_custom_y[il][2*tuple] - exact_sol_y[il][2*tuple]),
                                                std::abs(bin_custom_y[il][2*tuple] - exact_sol_y[il][2*tuple])/exact_sol_y[il][2*tuple]);
        if( dist > tol )
        {
         std::cout << std::scientific << std::setprecision(16)
                   << "Error: Mismatch in bin values for " << testname                                                      << std::endl
                   << "case ("            << bin_custom_x[il][2*tuple]   << ";"   << bin_custom_x[il + 1][2*tuple]   << ")" << std::endl 
                   << "bin = "            << bin_custom_y[il][2*tuple]                                                      << std::endl
                   << "bin_exact = "      << exact_sol_y[il][2*tuple]                                                       << std::endl
                   << "relative error = " << dist                                                                           << std::endl
                   << "tolerance = "      << tol                                                                            << std::endl;

         return_flag = 1;
        }

//tuple + 1
        dist = Antioch::if_else(exact_sol_y[il][2*tuple + 1] < tol,
                                                std::abs(bin_custom_y[il][2*tuple + 1] - exact_sol_y[il][2*tuple + 1]),
                                                std::abs(bin_custom_y[il][2*tuple + 1] - exact_sol_y[il][2*tuple + 1])/exact_sol_y[il][2*tuple + 1]);
        if( dist > tol )
        {
         std::cout << std::scientific << std::setprecision(16)
                   << "Error: Mismatch in bin values for " << testname                                                      << std::endl
                   << "case ("            << bin_custom_x[il][2*tuple + 1]   << ";"   << bin_custom_x[il+1][2*tuple + 1]   << ")" << std::endl 
                   << "bin = "            << bin_custom_y[il][2*tuple + 1]                                                        << std::endl
                   << "bin_exact = "      << exact_sol_y[il][2*tuple + 1]                                                         << std::endl
                   << "relative error = " << dist                                                                                 << std::endl
                   << "tolerance = "      << tol                                                                                  << std::endl;

         return_flag = 1;
        }


      }
    }
  }

  return return_flag;
}

int main()
{

  int returnval(0);

  returnval = returnval ||
    vectester (std::valarray<float>(2*ANTIOCH_N_TUPLES), "valarray<float>");
  returnval = returnval ||
    vectester (std::valarray<double>(2*ANTIOCH_N_TUPLES), "valarray<double>");
  returnval = returnval ||
    vectester (std::valarray<long double>(2*ANTIOCH_N_TUPLES), "valarray<ld>");
/*
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
*/
#ifdef ANTIOCH_HAVE_GRVY
  gt.Finalize();
  gt.Summarize();
#endif

 return returnval;

}
