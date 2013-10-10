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
  std::ifstream hv(path_to_files + "/solar_flux.dat");

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
  }
  CH4.close();

  while(!hv.eof())
  {
    Scalar w,l,dw;
    hv >> l >> w >> dw;
    hv_lambda.push_back(l*10.L); //nm -> Angström
    hv_irr.push_back(w * 1e4 * Antioch::Constants::Planck_constant<Scalar>() 
                             * Antioch::Constants::light_celerity<Scalar>());//m-2 -> cm-2 and times h*c/lamdba energy -> number of photons
  }
  hv.close();

  Scalar T(1500.);
  Antioch::ParticleFlux<std::vector<Scalar> > hv_flux(hv_lambda,hv_irr);

  Antioch::PhotochemicalRate<Scalar, std::vector<Scalar> > rate_hv(CH4_cs,CH4_lambda);

  Scalar rate_exact = 2.5838962773224599416e-40;
  Antioch::SigmaBinConverter<std::vector<Scalar> > youpi;
  std::vector<Scalar> irr_rescaled,sigma_rescaled;

  rate_hv.calculate_rate_constant(hv_flux.flux(), hv_flux.abscissa());
  Scalar rate = rate_hv.rate(T);
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;

  int return_flag = 0;

  if( abs( (rate - rate_exact)/rate_exact ) > tol )
  {
    std::cout << std::scientific << std::setprecision(16)
              << "Error: Mismatch in rate values." << std::endl
              << "rate = "       << rate           << std::endl
              << "rate_exact = " << rate_exact     << std::endl;

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

  return (tester<float>(std::string(argv[1])) ||
          tester<double>(std::string(argv[1])) ||
          tester<long double>(std::string(argv[1]))
          );
}
