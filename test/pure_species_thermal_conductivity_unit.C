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
#include <cmath>

// Antioch
#include "antioch/stat_mech_thermo.h"
#include "antioch/pure_species_thermal_conductivity.h"
#include "antioch/chemical_mixture.h"
#include "antioch/metaprogramming.h"

template <typename Scalar>
int test_k( const Scalar k, const Scalar k_exact, const Scalar tol )
{
  using std::abs;

  int return_flag = 0;

  const Scalar rel_error = abs( (k - k_exact)/k_exact);

  if( rel_error  > tol )
    {
      std::cerr << std::setprecision(15)
                << "Error: Mismatch in thermal conductivity" << std::endl
		<< "k       = " << k << std::endl
		<< "k_exact = " << k_exact << std::endl
		<< "rel_error = " << rel_error << std::endl
		<< "tol = " << tol << std::endl;
      return_flag = 1;
    }

  return return_flag;
}

template <typename Scalar>
int tester()
{
  std::vector<std::string> species_str_list;
  const unsigned int n_species = 1;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );

  const Scalar LJ_depth_N2 = 97.53L;
  const Scalar Z_298 = 4.0L;

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );

  Antioch::StatMechThermodynamics<Scalar> thermo( chem_mixture );

  Antioch::PureSpeciesThermalConductivity<Antioch::StatMechThermodynamics<Scalar>, Scalar > k( thermo, Z_298, LJ_depth_N2);

  const Scalar mu  = 3.14e-3;
  const Scalar dss = 5.23e-5;
  const Scalar rho = 1.4;
  const Scalar T = 1500.1L;

  // from bc
  const Scalar k_N2_exact = 3.194342919259960202334421163642706718735099613817392359646;
  int return_flag = 0;

  const Scalar tol = (std::numeric_limits<Scalar>::epsilon() * 2 < 7e-17)?7e-17:
                                                                           std::numeric_limits<Scalar>::epsilon() * 2;


  int return_flag_temp = 0;
  return_flag_temp = test_k( k(0,mu,T,rho,dss), k_N2_exact, tol );
  if( return_flag_temp != 0 ) return_flag = 1;

  return return_flag;
}

int main()
{
  return (tester<double>() ||
         tester<long double>() ||
          tester<float>());
}
