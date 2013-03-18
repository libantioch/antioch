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
#include <cmath>
#include <limits>

// Antioch
#include "antioch/physical_constants.h"
#include "antioch/chemical_mixture.h"
#include "antioch/stat_mech_thermo.h"

template <typename Scalar>
bool test_relative(const Scalar val, const Scalar truth, const Scalar tol)
{
  if( std::abs( (val-truth)/truth ) > tol )
    return false;
  else
    return true;
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

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );

  // Can we instantiate it?
  Antioch::StatMechThermodynamics<Scalar> sm_thermo( chem_mixture );

  // //const Scalar P = 100000.0;
  // const std::vector<Scalar> mass_fractions( 5, 0.2 );
  // const Scalar T1 = 190.0;
  // const Scalar T2 = 1500.0;
  // const Scalar T3 = 10000.0;

  const Scalar R_N2 = Antioch::Constants::R_universal/28.016;
  const Scalar R_O2 = Antioch::Constants::R_universal/32.0;
  const Scalar R_N = Antioch::Constants::R_universal/14.008;
  const Scalar R_O = Antioch::Constants::R_universal/16.0;
  const Scalar R_NO = Antioch::Constants::R_universal/30.008;

  int return_flag = 0;

  Scalar tol = std::numeric_limits<Scalar>::epsilon() * 2;

  // N2
  {
    Scalar cv_N2 = sm_thermo.cv_tr(0);
    
    if( !test_relative(cv_N2, R_N2*Scalar(2.5), tol) )
    {
      std::cerr << std::scientific << std::setprecision(20);
      std::cerr << "Error: Mismatch in cv_tr for N2."
                << "\n Expected = " << R_N2*Scalar(2.5)
                << "\n Computed = " << cv_N2
                << "\n Diff     = " << cv_N2 - R_N2*Scalar(2.5)
                << std::endl;
      return_flag += 1;
    }
  }

  // O2
  {
    Scalar cv_O2 = sm_thermo.cv_tr(1);
    
    if( !test_relative(cv_O2, R_O2*Scalar(2.5), tol) )
    {
      std::cerr << std::scientific << std::setprecision(20);
      std::cerr << "Error: Mismatch in cv_tr for O2."
                << "\n Expected = " << R_O2*Scalar(2.5)
                << "\n Computed = " << cv_O2
                << "\n Diff     = " << cv_O2 - R_O2*Scalar(2.5)
                << std::endl;
      return_flag += 1;
    }
  }

  // N
  {
    Scalar cv_N = sm_thermo.cv_tr(2);
    
    if( !test_relative(cv_N, R_N*Scalar(1.5), tol) )
    {
      std::cerr << std::scientific << std::setprecision(20);
      std::cerr << "Error: Mismatch in cv_tr for N."
                << "\n Expected = " << R_N*Scalar(1.5)
                << "\n Computed = " << cv_N
                << "\n Diff     = " << cv_N - R_N*Scalar(2.5)
                << std::endl;
      return_flag += 1;
    }
  }

  // O
  {
    Scalar cv_O = sm_thermo.cv_tr(3);
    
    if( !test_relative(cv_O, R_O*Scalar(1.5), tol) )
    {
      std::cerr << std::scientific << std::setprecision(20);
      std::cerr << "Error: Mismatch in cv_tr for O."
                << "\n Expected = " << R_O*Scalar(1.5)
                << "\n Computed = " << cv_O
                << "\n Diff     = " << cv_O - R_O*Scalar(2.5)
                << std::endl;
      return_flag += 1;
    }
  }

  // NO
  {
    Scalar cv_NO = sm_thermo.cv_tr(4);
    
    if( !test_relative(cv_NO, R_NO*Scalar(2.5), tol) )
    {
      std::cerr << std::scientific << std::setprecision(20);
      std::cerr << "Error: Mismatch in cv_tr for NO."
                << "\n Expected = " << R_NO*Scalar(2.5)
                << "\n Computed = " << cv_NO
                << "\n Diff     = " << cv_NO - R_NO*Scalar(2.5)
                << std::endl;
      return_flag += 1;
    }
  }

  return return_flag;
}


int main()
{
// We're not getting the full long double precision yet?
  return (tester<double>() ||
//          tester<long double>() ||
          tester<float>());
}
