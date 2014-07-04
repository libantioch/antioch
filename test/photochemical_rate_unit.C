//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014 Paul T. Bauman, Sylvain Plessis, Roy H. Stonger
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

#include "antioch/particle_flux.h"
#include "antioch/photochemical_rate.h"
#include "antioch/physical_constants.h"

#include "antioch/vector_utils.h"

template <typename Scalar>
int tester(std::string path_to_files)
{
  std::ifstream CH4(path_to_files + "/CH4_hv_cs.dat");
  std::ifstream  hv(path_to_files + "/solar_flux.dat");

  std::string first_line;

  getline(CH4,first_line);
  getline(hv,first_line);

  std::vector<Scalar> CH4_cs;
  std::vector<Scalar> CH4_lambda;
  std::vector<Scalar> hv_irr;
  std::vector<Scalar> hv_lambda;

  while(!CH4.eof())
  {
    Scalar cs,l;
    CH4 >> l >> cs;
    CH4_lambda.push_back(l);
    CH4_cs.push_back(cs);
    if(CH4_lambda.size() == 137)break;
  }
  CH4.close();

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

  Scalar T(1500.L);

  Antioch::PhotochemicalRate<Scalar, std::vector<Scalar> > rate_hv(CH4_cs,CH4_lambda);

  Antioch::SigmaBinConverter<std::vector<Scalar> > bin;
  std::vector<Scalar> sigma_rescaled;
  bin.y_on_custom_grid(CH4_lambda,CH4_cs,hv_lambda,sigma_rescaled);

  rate_hv.calculate_rate_constant(hv_irr, hv_lambda);
  Scalar rate = rate_hv.rate(T);

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;
  Scalar rate_exact(0.L);

  for(unsigned int il = 0; il < hv_lambda.size() - 1; il++)
  {
      rate_exact += sigma_rescaled[il] * hv_irr[il] * (hv_lambda[il+1] - hv_lambda[il]);
  }

  int return_flag = (rate_exact == Scalar(0.L));
  if(return_flag)std::cout << "Error: rate is null" << std::endl;
  if( std::abs( (rate - rate_exact)/rate_exact ) > tol )
  {
    std::cout << std::scientific << std::setprecision(16)
              << "Error: Mismatch in rate values."     << std::endl
              << "rate = "           << rate           << std::endl
              << "rate_exact = "     << rate_exact     << std::endl
              << "relative error = " << std::abs(rate_exact - rate)/rate_exact << std::endl
              << "tolerance = "      << tol            << std::endl;

    return_flag = 1;
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
