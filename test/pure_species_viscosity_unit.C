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

// C++
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

// Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/physical_constants.h"
#include "antioch/pure_species_viscosity.h"
#include "antioch/vector_utils.h"

template <typename Scalar>
int test_viscosity( const Scalar mu, const Scalar mu_exact, const Scalar tol )
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

template <typename Scalar>
int tester()
{
// value for N2
  const Scalar LJ_depth(97.530L);
  const Scalar LJ_diameter(3.621L);
  const Scalar dipole_moment(0.L);
  const Scalar mass(28.016e-3L/Antioch::Constants::Avogadro<Scalar>());

  Antioch::PureSpeciesViscosity<Scalar,Antioch::GSLSpliner> mu(LJ_depth,LJ_diameter,dipole_moment,mass);

  const Scalar T = 1500.1;

  // bc gives
  const Scalar mu_exact_times_interp = 0.0000417395098853601937871105407365424874568203066945573066;

  int return_flag = 0;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 10.;

  return_flag = test_viscosity( mu(T) * mu.Stockmayer(T), mu_exact_times_interp, tol );

  return return_flag;
}

int main()
{
  return tester<double>() ||
         tester<long double>() ||
         tester<float>();
};
