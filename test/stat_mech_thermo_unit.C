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
bool test_zero(const Scalar val, const Scalar tol)
{
  if( std::abs(val) > tol )
    return false;
  else
    return true;
}

template <typename Scalar>
int test_cv_tr()
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

  // Mass fractions
  std::vector<Scalar> mass_fractions( 5, 0.2 );
  mass_fractions[0] = 0.5;
  mass_fractions[1] = 0.2;
  mass_fractions[2] = 0.1;
  mass_fractions[3] = 0.1;
  mass_fractions[4] = 0.1;

  Scalar cv_tr_mix = 0.0;

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

    cv_tr_mix += cv_N2*mass_fractions[0];
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

    cv_tr_mix += cv_O2*mass_fractions[1];
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

    cv_tr_mix += cv_N*mass_fractions[2];
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
    
    cv_tr_mix += cv_O*mass_fractions[3];
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

    cv_tr_mix += cv_NO*mass_fractions[4];
  }

  // mixture
  {
    Scalar cv = sm_thermo.cv_tr(mass_fractions);
    
    if( !test_relative(cv, cv_tr_mix, tol) )
    {
      std::cerr << std::scientific << std::setprecision(20);
      std::cerr << "Error: Mismatch in mixture cv_tr."
                << "\n Expected = " << cv_tr_mix
                << "\n Computed = " << cv
                << "\n Diff     = " << cv - cv_tr_mix
                << std::endl;
      return_flag += 1;
    }
  }

  return return_flag;
}

template <typename Scalar>
int test_cv_vib()
{
  using std::exp;

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

  // Mass fractions
  std::vector<Scalar> mass_fractions( 5, 0.2 );
  mass_fractions[0] = 0.5;
  mass_fractions[1] = 0.2;
  mass_fractions[2] = 0.1;
  mass_fractions[3] = 0.1;
  mass_fractions[4] = 0.1;

  const Scalar R_N2 = Antioch::Constants::R_universal/28.016;
  const Scalar R_O2 = Antioch::Constants::R_universal/32.0;
  const Scalar R_N = Antioch::Constants::R_universal/14.008;
  const Scalar R_O = Antioch::Constants::R_universal/16.0;
  const Scalar R_NO = Antioch::Constants::R_universal/30.008;

  const Scalar th0_N2 = 3.39500e+03; // degeneracy = 1
  const Scalar th0_O2 = 2.23900e+03; // degeneracy = 1
  const Scalar th0_NO = 2.81700e+03; // degeneracy = 1

  // Tv
  const Scalar Tv = 1000.0;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 2;
  const Scalar ztol = std::numeric_limits<Scalar>::epsilon();

  int return_flag = 0;

  Scalar cv_vib_mix = 0.0;

  // N2
  {
    Scalar cv_vib_N2 = sm_thermo.cv_vib (0, Tv);

    const Scalar expv   = exp(th0_N2/Tv);
    const Scalar expvmi = expv - Scalar(1.0);
    Scalar cv_vib_N2_true = R_N2*th0_N2*th0_N2*expv/expvmi/expvmi/Tv/Tv;

    if( !test_relative(cv_vib_N2, cv_vib_N2_true, tol) )
    {
      std::cerr << std::scientific << std::setprecision(20);
      std::cerr << "Error: Mismatch in cv_vib for N2."
                << "\n Expected = " << cv_vib_N2_true
                << "\n Computed = " << cv_vib_N2
                << "\n Diff     = " << cv_vib_N2_true - cv_vib_N2
                << std::endl;
      return_flag += 1;
    }

    cv_vib_mix += mass_fractions[0]*cv_vib_N2_true;
  }

  // O2
  {
    Scalar cv_vib_O2 = sm_thermo.cv_vib (1, Tv);

    const Scalar expv   = exp(th0_O2/Tv);
    const Scalar expvmi = expv - Scalar(1.0);
    Scalar cv_vib_O2_true = R_O2*th0_O2*th0_O2*expv/expvmi/expvmi/Tv/Tv;

    if( !test_relative(cv_vib_O2, cv_vib_O2_true, tol) )
    {
      std::cerr << std::scientific << std::setprecision(20);
      std::cerr << "Error: Mismatch in cv_vib for O2."
                << "\n Expected = " << cv_vib_O2_true
                << "\n Computed = " << cv_vib_O2
                << "\n Diff     = " << cv_vib_O2_true - cv_vib_O2
                << std::endl;
      return_flag += 1;
    }

    cv_vib_mix += mass_fractions[1]*cv_vib_O2_true;
  }

  // O
  {
    Scalar cv_vib_O = sm_thermo.cv_vib (2, Tv);

    if( !test_zero(cv_vib_O, ztol) )
    {
      std::cerr << std::scientific << std::setprecision(20);
      std::cerr << "Error: Mismatch in cv_vib for O."
                << "\n Expected = " << Scalar(0.0)
                << "\n Computed = " << cv_vib_O
                << "\n Diff     = " << cv_vib_O
                << std::endl;
      return_flag += 1;
    }

    // cv_vib_mix += 0.0;
  }

  // N
  {
    Scalar cv_vib_N = sm_thermo.cv_vib (3, Tv);

    if( !test_zero(cv_vib_N, ztol) )
    {
      std::cerr << std::scientific << std::setprecision(20);
      std::cerr << "Error: Mismatch in cv_vib for N."
                << "\n Expected = " << Scalar(0.0)
                << "\n Computed = " << cv_vib_N
                << "\n Diff     = " << cv_vib_N
                << std::endl;
      return_flag += 1;
    }

    // cv_vib_mix += 0.0;
  }

  // NO
  {
    Scalar cv_vib_NO = sm_thermo.cv_vib (4, Tv);

    const Scalar expv   = exp(th0_NO/Tv);
    const Scalar expvmi = expv - Scalar(1.0);
    Scalar cv_vib_NO_true = R_NO*th0_NO*th0_NO*expv/expvmi/expvmi/Tv/Tv;

    if( !test_relative(cv_vib_NO, cv_vib_NO_true, tol) )
    {
      std::cerr << std::scientific << std::setprecision(20);
      std::cerr << "Error: Mismatch in cv_vib for NO."
                << "\n Expected = " << cv_vib_NO_true
                << "\n Computed = " << cv_vib_NO
                << "\n Diff     = " << cv_vib_NO_true - cv_vib_NO
                << std::endl;
      return_flag += 1;
    }

    cv_vib_mix += mass_fractions[4]*cv_vib_NO_true;
  }

  // mixture
  {
    Scalar cv = sm_thermo.cv_vib(Tv, mass_fractions);
    
    if( !test_relative(cv, cv_vib_mix, tol) )
    {
      std::cerr << std::scientific << std::setprecision(20);
      std::cerr << "Error: Mismatch in mixture cv_vib."
                << "\n Expected = " << cv_vib_mix
                << "\n Computed = " << cv
                << "\n Diff     = " << cv - cv_vib_mix
                << std::endl;
      return_flag += 1;
    }
  }

  return return_flag;
}

template <typename Scalar>
int test_cv_el()
{
  using std::exp;

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

  // Mass fractions
  std::vector<Scalar> mass_fractions( 5, 0.2 );
  mass_fractions[0] = 0.5;
  mass_fractions[1] = 0.2;
  mass_fractions[2] = 0.1;
  mass_fractions[3] = 0.1;
  mass_fractions[4] = 0.1;

  const Scalar R_N2 = Antioch::Constants::R_universal/28.016;
  const Scalar R_O2 = Antioch::Constants::R_universal/32.0;
  const Scalar R_N = Antioch::Constants::R_universal/14.008;
  const Scalar R_O = Antioch::Constants::R_universal/16.0;
  const Scalar R_NO = Antioch::Constants::R_universal/30.008;

  // Te
  const Scalar Te = 1000.0;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 2;
  const Scalar ztol = std::numeric_limits<Scalar>::epsilon();

  int return_flag = 0;

  // N2
  {
    Scalar cv_el_N2 = sm_thermo.cv_el (0, Te);
  }

  // Mixture
  {
    Scalar cv_el = sm_thermo.cv_el (Te, mass_fractions);
  }

  // TODO: Add some non-trivial checks

  return return_flag;
}

int main()
{

  // Translational/rotational specific heat at constant volume
  int ierr = (test_cv_tr<double>() ||
//            test_cv_tr<long double>() ||
              test_cv_tr<float>());

  // Vibrational specific heat at constant volume
  ierr += (test_cv_vib<double>() ||
//         test_cv_vib<long double>() ||
           test_cv_vib<float>());

  // Electronic specific heat at constant volume
  ierr += (test_cv_el<double>() ||
//         test_cv_el<long double>() ||
           test_cv_el<float>());

  return ierr;
}
