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

// Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/gsl_spliner.h"
#include "antioch/vector_utils.h"

// C++
#include <limits>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>

template <typename Scalar>
int check_value(const Scalar & ref, const Scalar & candidate, const Scalar & x, const std::string & words)
{
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 10.;

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
  return 10.L + 5.L * x + 10.L * x * x - 2.L * x * x * x;
}

template <typename Scalar>
void fill_ref(std::vector<Scalar> & x_ref, std::vector<Scalar> & y_ref,
              unsigned int n_data, const Scalar & min,
              const Scalar & max)
{
std::ofstream out("gsl.ref");
  for(unsigned int i = 0; i < n_data; i++)
  {
     x_ref[i] = min + (Scalar)(i) * (max - min) / (Scalar)(n_data-1);
     y_ref[i] = function(x_ref[i]);
out << x_ref[i] << " " << y_ref[i] << std::endl;
  }
  out.close();

}

template <typename Scalar>
int tester()
{
  const unsigned int n_data(40);
  const unsigned int n_test(15);
  std::vector<Scalar> x_ref(n_data,0),y_ref(n_data,0);

  const Scalar min = -5L;
  const Scalar max = 8L;
  fill_ref(x_ref,y_ref,n_data,min, max);

  Antioch::GSLSpliner default_constructor;
  Antioch::GSLSpliner explicit_constructor(x_ref,y_ref);

  default_constructor.spline_init(x_ref,y_ref);

  gsl_interp_accel * acc;
  gsl_spline       * spline;
  acc =  gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline, x_ref.size());
  double xptr[n_data]; 
  double yptr[n_data]; 
  for(unsigned int i = 0; i < n_data; ++i)
  {
      xptr[i] = x_ref[i];
      yptr[i] = y_ref[i];
  }
  gsl_spline_init(spline, xptr, yptr, n_data);

  int return_flag(0);

  std::ofstream out("gsl.dat");
  for(unsigned int n = 0; n < n_test; n++)
  {
     Scalar x = min + (Scalar)(n) * (max - min) / (Scalar)(n_test-1);
     Scalar exact = function(x);
     Scalar spline_default = default_constructor.interpolated_value(x);
     Scalar spline_explicit = explicit_constructor.interpolated_value(x);
     Scalar  y = gsl_spline_eval(spline,x,acc);
        out << x << " " << exact << " " << y << std::endl;
//     return_flag = check_value(exact,spline_default,"default constructor") || return_flag;
 //    return_flag = check_value(exact,spline_explicit,"explicit constructor") || return_flag;
     return_flag = check_value(exact,spline_default,x,"GODDAMMIT") || return_flag;
  }
  out.close();


  return return_flag;
}

int main()
{
// gsl work in double...
  return (tester<float>() ||
          tester<double>());
 //         tester<long double>() ||
}
