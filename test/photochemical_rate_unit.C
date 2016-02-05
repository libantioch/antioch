//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014-2016 Paul T. Bauman, Benjamin S. Kirk,
//                         Sylvain Plessis, Roy H. Stonger
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
#include "antioch/kinetics_parsing.h"

#include "antioch/vector_utils.h"

template <typename Scalar>
int check_rate(const Scalar & rate_exact, const Scalar & rate)
{
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;
  int return_flag = (rate_exact > tol); // not zero
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
    Scalar cs,l(-1);
    CH4 >> l >> cs;
    if(!CH4.good())break;
    CH4_lambda.push_back(l);
    CH4_cs.push_back(cs);
  }
  CH4.close();

  while(!hv.eof())
  {
    Scalar w,l(-1),dw;
    hv >> l >> w >> dw;
    if(!hv.good())break;
    hv_lambda.push_back(l * 10); //nm -> Angström
    hv_irr.push_back(w * 1e-4L  // * 1e-4: m-2 -> cm-2 
                       / (Antioch::Constants::Planck_constant<Scalar>() * Antioch::Constants::light_celerity<Scalar>() / l)// /(h*c/lambda): energy -> number of photons
                       / 10); // by Angström
  }
  hv.close();

  Antioch::ParticleFlux<std::vector<Scalar> > part_flux(hv_lambda,hv_irr);

  Antioch::PhotochemicalRate<Scalar, std::vector<Scalar> > rate_hv(CH4_cs,CH4_lambda);

  Antioch::SigmaBinConverter<std::vector<Scalar> > bin;
  std::vector<Scalar> sigma_rescaled(hv_lambda.size());
  bin.y_on_custom_grid(CH4_lambda,CH4_cs,hv_lambda,sigma_rescaled);

  Scalar rate = rate_hv.rate(part_flux);

  Scalar rate_exact(0);

  for(unsigned int il = 0; il < hv_lambda.size() - 1; il++)
  {
      rate_exact += sigma_rescaled[il] * hv_irr[il] * (hv_lambda[il+1] - hv_lambda[il]);
  }

  int return_flag = check_rate(rate_exact,rate);

 // multiplying by 2 the cross-section
  int il = CH4_cs.size() * 2 / 3; 
  CH4_cs[il] *= 2;

  bin.y_on_custom_grid(CH4_lambda,CH4_cs,hv_lambda,sigma_rescaled);

  Antioch::set_zero(rate_exact);
  for(unsigned int il = 0; il < hv_lambda.size() - 1; il++)
  {
      rate_exact += sigma_rescaled[il] * hv_irr[il] * (hv_lambda[il+1] - hv_lambda[il]);
  }
  rate_hv.set_parameter(Antioch::KineticsModel::Parameters::SIGMA, il, CH4_cs[il]);
  rate = rate_hv.rate(part_flux);

  return_flag = check_rate(rate_exact,rate) || return_flag;

 // multiplying by 2 one value of the cross-section
  il = CH4_cs.size()/2;
  CH4_cs[il] *= 2;

  bin.y_on_custom_grid(CH4_lambda,CH4_cs,hv_lambda,sigma_rescaled);

  Antioch::set_zero(rate_exact);
  for(unsigned int il = 0; il < hv_lambda.size() - 1; il++)
  {
      rate_exact += sigma_rescaled[il] * hv_irr[il] * (hv_lambda[il+1] - hv_lambda[il]);
  }

// the other way to reset
  Antioch::reset_parameter_of_rate(rate_hv,Antioch::KineticsModel::Parameters::SIGMA, CH4_cs[il] , il, "SI");
  rate = rate_hv.rate(part_flux);

  return_flag = check_rate(rate_exact,rate) || return_flag;

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
