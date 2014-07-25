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

template <typename Scalar>
int tester(std::string path_to_files)
{
  std::ifstream  hv(path_to_files + "/solar_flux.dat");

  std::string first_line;

  getline(hv,first_line);

  std::vector<Scalar> hv_irr;
  std::vector<Scalar> hv_lambda;

  while(!hv.eof())
  {
    Scalar w,l,dw;
    hv >> l >> w >> dw;
    hv_lambda.push_back(l * 10.L); //nm -> Angström
    hv_irr.push_back(w * 1e-4L  // * 1e-4: m-2 -> cm-2 
                       / (Antioch::Constants::Planck_constant<Scalar>() * Antioch::Constants::light_celerity<Scalar>() / l)// /(h*c/lambda): energy -> number of photons
                       / 10.); // by Angström
    if(hv_lambda.size() == 796)break;
  }
  hv.close();

  Antioch::SigmaBinConverter<std::vector<Scalar> > bin;
  std::vector<Scalar> flux_recalculated;
  bin.y_on_custom_grid(hv_lambda,hv_irr,hv_lambda,flux_recalculated);

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;

  int return_flag = 0;
  for(unsigned int il = 0; il < hv_lambda.size() - 1; il++)
  {
    if( std::abs( (flux_recalculated[il] - hv_irr[il])/hv_irr[il]) > tol )
    {
      std::cout << std::scientific << std::setprecision(16)
                << "Error: Mismatch in bin values."            << std::endl
                << "bin = "           << flux_recalculated[il] << std::endl
                << "bin_exact = "     << hv_irr[il]            << std::endl
                << "relative error = " << std::abs( (flux_recalculated[il] - hv_irr[il])/hv_irr[il]) << std::endl
                << "tolerance = "      << tol            << std::endl;

      return_flag = 1;
    }
  }

  return return_flag;
}

int main(int argc, char* argv[])
{
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify input files location." << std::endl;
      antioch_error();
    }

   return (tester<float>(std::string(argv[1]))  ||
           tester<double>(std::string(argv[1])) ||
           tester<long double>(std::string(argv[1]))
          );
}
