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

// C++
#include <limits>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

// Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/sigma_bin_converter.h"
#include "antioch/antioch_asserts.h"
#include "antioch/physical_constants.h"
#include "antioch/vector_utils.h"
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
template <typename VectorScalar>
void make_custom(unsigned int choice, VectorScalar & x, VectorScalar & y)
{
  // sum_{bin} Delta x * y
  switch(choice)
  {
      case(0): // 9 bin contained in ref -> [2.5;8.5] within [1;10]
      {
          x.resize(9);   y.resize(9);
          x[0] = 2.50L;  y[0] = 2.25L/0.75L; // 0.50 * 2  + 0.25 * 5
          x[1] = 3.25L;  y[1] = 3.75L/0.75L; // 0.75 * 5
          x[2] = 4.00L;  y[2] = 6.00L/0.75L; // 0.75 * 8
          x[3] = 4.75L;  y[3] = 5.00L/0.75L; // 0.25 * 8  + 0.5  * 6
          x[4] = 5.50L;  y[4] = 5.50L/0.75L; // 0.25 * 6  + 0.25 * 10
          x[5] = 6.25L;  y[5] = 7.50L/0.75L; // 0.75 * 10
          x[6] = 7.00L;  y[6] = 5.25L/0.75L; // 0.75 * 7
          x[7] = 7.75L;  y[7] = 3.75L/0.75L; // 0.25 * 7 + 0.5 * 4
          x[8] = 8.50L;  y[8] = 0.00L/0.75L; // 0. (right stairs)
        break;
      }
      case(1):// 9 bin outside ref -> [0;12] containing [1;10]
      {
          x.resize(9);   y.resize(9);
          x[0] = 0.00L;  y[0] = 0.50L/1.50L; // 0.5 * 1
          x[1] = 1.50L;  y[1] = 2.50L/1.50L; // 0.5 * 1  + 2.0 * 1
          x[2] = 3.00L;  y[2] = 9.00L/1.50L; // 1.0 * 5  + 0.5 * 8
          x[3] = 4.50L;  y[3] = 10.0L/1.50L; // 0.5 * 8  + 1.0 * 6
          x[4] = 6.00L;  y[4] = 13.5L/1.50L; // 1.0 * 10 + 0.5 * 7
          x[5] = 7.50L;  y[5] = 7.50L/1.50L; // 0.50 * 7 + 1.0 * 4
          x[6] = 9.00L;  y[6] = 2.00L/1.50L; // 1.0 * 2 
          x[7] = 10.5L;  y[7] = 0.00L/1.50L; // 0
          x[8] = 12.0L;  y[8] = 0.00L/1.50L; // 0. (right stairs)
        break;
      }
      case(2): // 5 bins beyond only min -> [-1;3.8] in [0;10]
      {
          x.resize(5);    y.resize(5);
          x[0] = -1.00L;  y[0] = 0.00L/1.20L; // 0
          x[1] =  0.20L;  y[1] = 0.40L/1.20L; // 0.4 * 1
          x[2] =  1.40L;  y[2] = 1.80L/1.20L; // 0.6 * 1  + 0.6 * 2
          x[3] =  2.60L;  y[3] = 4.80L/1.20L; // 0.4 * 2  + 0.8 * 5
          x[4] =  3.80L;  y[4] = 00.0L/1.20L; // 0. (right stairs)
        break;
      }
      case(3):// 6 bins beyond only max -> [2;10.75] in [0;10]
      {
          x.resize(6);    y.resize(6);
          x[0] =  2.00L;  y[0] =  5.75L/1.75L; // 1.00 * 2 + 0.75 * 5
          x[1] =  3.75L;  y[1] = 12.25L/1.75L; // 0.25 * 5 + 1.00 * 8 + 0.50 * 6
          x[2] =  5.50L;  y[2] = 14.75L/1.75L; // 0.50 * 6 + 1.00 * 10 + 0.25 * 7
          x[3] =  7.25L;  y[3] =  9.25L/1.75L; // 0.75 * 7 + 1.00 * 4
          x[4] =  9.00L;  y[4] =  2.00L/1.75L; // 1.00 * 2
          x[5] = 10.75L;  y[5] = 00.00L/1.75L; // 0. (right stairs)
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


template <typename Scalar>
int tester()
{
  std::vector<Scalar> bin_ref_x,bin_ref_y;

  make_reference(bin_ref_x,bin_ref_y);

  Antioch::SigmaBinConverter<std::vector<Scalar> > bin;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 10;

  int return_flag = 0;

  // 4 cases:
  //   - custom inside ref
  //   - ref inside custom
  //   - custom beyond only min ref
  //   - custom beyond only max ref
  for(unsigned int i = 0; i < 4; i++)
  {
    std::vector<Scalar> bin_custom_x, exact_sol_y;
    std::vector<Scalar> bin_custom_y;
    make_custom(i,bin_custom_x,exact_sol_y);

    bin.y_on_custom_grid(bin_ref_x,bin_ref_y,
                         bin_custom_x,bin_custom_y);


    for(unsigned int il = 0; il < bin_custom_x.size() - 1; il++)
    {
      const Scalar dist = (exact_sol_y[il] < tol)?std::abs(bin_custom_y[il] - exact_sol_y[il]):std::abs(bin_custom_y[il] - exact_sol_y[il])/exact_sol_y[il];
      if( dist > tol )
      {
        std::cout << std::scientific << std::setprecision(16)
                  << "Error: Mismatch in bin values."                                                  << std::endl
                  << "case ("            << bin_custom_x[il]   << ";"   << bin_custom_x[il+1]   << ")" << std::endl 
                  << "bin = "            << bin_custom_y[il]                                           << std::endl
                  << "bin_exact = "      << exact_sol_y[il]                                            << std::endl
                  << "relative error = " << dist                                                       << std::endl
                  << "tolerance = "      << tol                                                        << std::endl;

        return_flag = 1;
      }
    }
  }

  return return_flag;
}

int main()
{

   return (tester<float>()  ||
           tester<double>() ||
           tester<long double>()
          );
}
