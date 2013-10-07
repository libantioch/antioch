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
#include <iostream>
#include <cmath>

// Antioch
#include "antioch/stat_mech_thermo.h"
#include "antioch/eucken_thermal_conductivity.h"

template <typename Scalar>
int test_k( const Scalar k, const Scalar k_exact, const Scalar tol )
{
  using std::abs;

  int return_flag = 0;

  const Scalar rel_error = abs( (k - k_exact)/k_exact);

  if( rel_error  > tol )
    {
      std::cerr << "Error: Mismatch in thermal conductivity" << std::endl
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
  const unsigned int n_species = 5;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "O2" );
  species_str_list.push_back( "N" );
  species_str_list.push_back( "O" );
  species_str_list.push_back( "NO" );

  const Scalar Mm_N = 14.008e-3L;
  const Scalar Mm_O = 16.000e-3L;
  const Scalar Mm_N2 = 2.L * Mm_N;
  const Scalar Mm_O2 = 2.L * Mm_O;
  const Scalar Mm_NO = Mm_N + Mm_O;

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );

  const Scalar R_N2 = Antioch::Constants::R_universal<Scalar>()/Mm_N2;
  const Scalar R_O2 = Antioch::Constants::R_universal<Scalar>()/Mm_O2;
  const Scalar R_N = Antioch::Constants::R_universal<Scalar>()/Mm_N;
  const Scalar R_O = Antioch::Constants::R_universal<Scalar>()/Mm_O;
  const Scalar R_NO = Antioch::Constants::R_universal<Scalar>()/Mm_NO;

  Antioch::StatMechThermodynamics<Scalar> thermo( chem_mixture );

  Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<Scalar> > k( thermo );

  const Scalar mu = 3.14e-3;

  // octave gives
  const Scalar k_N2_trans_exact = 2.5*mu*1.5*R_N2;
  const Scalar k_N_trans_exact = 2.5*mu*1.5*R_N;
  const Scalar k_O2_trans_exact = 2.5*mu*1.5*R_O2;
  const Scalar k_O_trans_exact = 2.5*mu*1.5*R_O;
  const Scalar k_NO_trans_exact = 2.5*mu*1.5*R_NO;

  int return_flag = 0;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 2;

  int return_flag_temp = 0;
  return_flag_temp = test_k( k.trans(0, mu), k_N2_trans_exact, tol );
  if( return_flag_temp != 0 ) return_flag = 1;

  return_flag_temp = test_k( k.trans(1, mu), k_O2_trans_exact, tol );
  if( return_flag_temp != 0 ) return_flag = 1;

  return_flag_temp = test_k( k.trans(2, mu), k_N_trans_exact, tol );
  if( return_flag_temp != 0 ) return_flag = 1;

  return_flag_temp = test_k( k.trans(3, mu), k_O_trans_exact, tol );
  if( return_flag_temp != 0 ) return_flag = 1;

  return_flag_temp = test_k( k.trans(4, mu), k_NO_trans_exact, tol );
  if( return_flag_temp != 0 ) return_flag = 1;

  return return_flag;
}

int main()
{
  return (tester<double>() ||
//         tester<long double>() ||
          tester<float>());
}
