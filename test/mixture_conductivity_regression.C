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

// Antioch
#include "antioch_config.h"
#include "antioch/default_filename.h"
#include "antioch/vector_utils_decl.h"


#include "antioch/eucken_thermal_conductivity.h"
#include "antioch/eucken_thermal_conductivity_building.h"
#include "antioch/mixture_conductivity.h"

#include "antioch/sutherland_viscosity.h"
#include "antioch/blottner_viscosity.h"
#include "antioch/kinetics_theory_viscosity.h"
#include "antioch/mixture_viscosity.h"
#include "antioch/vector_utils.h"


#include "antioch/blottner_parsing.h"
#include "antioch/sutherland_parsing.h"
#include "antioch/kinetics_theory_viscosity_building.h"

#include "antioch/stat_mech_thermo.h"

template<typename Scalar>
bool test_val( Scalar value, Scalar gold, Scalar tol, std::string test )
{
  bool test_passed = true;

  Scalar error = std::abs( value - gold );
  if( error > tol )
    {
      test_passed = false;
      std::cerr << "Test "+test+" failed!" << std::endl
                << "value = " << value << std::endl
                << "gold  = " << gold << std::endl
                << "error = " << error << std::endl
                << "tol   = " << tol << std::endl;
    }

  return test_passed;
}

template <typename Scalar>
int tester()
{
  std::vector<std::string> species_str_list;

  // Test hard coded to 1 species right now. If we ever
  // add more here, be aware that the regression test on
  // the values is targeted at the N2 values, so if you
  // change the species order, update the N2_species index
  // accordingly.
  const unsigned int n_species = 1;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  unsigned int N2_species = 0;

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );

  typedef Antioch::StatMechThermodynamics<Scalar> MicroThermo;
  MicroThermo thermo_stat( chem_mixture );

  Antioch::TransportMixture<Scalar> tran_mixture( chem_mixture );

  Antioch::MixtureConductivity<Antioch::EuckenThermalConductivity<MicroThermo>,Scalar>
    e_k_mixture(tran_mixture);

  Antioch::build_eucken_thermal_conductivity<MicroThermo,Scalar>(e_k_mixture, thermo_stat);

  Antioch::MixtureViscosity<Antioch::SutherlandViscosity<Scalar>,Scalar>
    s_mu_mixture(tran_mixture);

  Antioch::MixtureViscosity<Antioch::BlottnerViscosity<Scalar>,Scalar>
    b_mu_mixture(tran_mixture);

#ifdef  ANTIOCH_HAVE_GSL
  Antioch::MixtureViscosity<Antioch::KineticsTheoryViscosity<Scalar, Antioch::GSLSpliner>,Scalar>
    k_mu_mixture(tran_mixture);
#endif // ANTIOCH_HAVE_GSL

  Antioch::read_sutherland_data_ascii<Scalar>( s_mu_mixture, Antioch::DefaultFilename::sutherland_data() );
  Antioch::read_blottner_data_ascii<Scalar>( b_mu_mixture, Antioch::DefaultFilename::blottner_data() );

#ifdef  ANTIOCH_HAVE_GSL
  Antioch::build_kinetics_theory_viscosity<Scalar>(k_mu_mixture);
#endif // ANTIOCH_HAVE_GSL

  const Scalar T = 1500.1;

  // Gold values were generated with long double version of test
  // Eucken, Sutherland gold
  Scalar k_e_s_N2_gold = 8.1294319021704618392660e-02L;
  Scalar k_e_s_N2_value = 0.0; // init

  // Eucken, Blottner gold
  Scalar k_e_b_N2_gold = 8.3906393746814049491975e-02L;
  Scalar k_e_b_N2_value = 0.0; // init

  // Eucken KineticsTheory gold
#ifdef  ANTIOCH_HAVE_GSL
  Scalar k_e_kt_N2_gold = 8.6774379182691310526782e-02L;
  Scalar k_e_kt_N2_value = 0.0; // init
#endif // ANTIOCH_HAVE_GSL

  std::cout << "Eucken (with SutherlandViscosity):" << std::endl;
  for( unsigned int s = 0; s < n_species; s++ )
    {
      Scalar mu = s_mu_mixture(s,T);
      Scalar value = e_k_mixture.conductivity_without_diffusion(s, T, mu);
      if( s == N2_species )
        k_e_s_N2_value = value;

      std::cout << "k(" << species_str_list[s] << ") = " << value << std::endl;
    }

  std::cout << "Eucken (with BlottnerViscosity):" << std::endl;
  for( unsigned int s = 0; s < n_species; s++ )
    {
      Scalar mu = b_mu_mixture(s,T);
      Scalar value = e_k_mixture.conductivity_without_diffusion(s, T, mu);
      if( s == N2_species )
        k_e_b_N2_value = value;

      std::cout << "k(" << species_str_list[s] << ") = " <<  value << std::endl;
    }

#ifdef  ANTIOCH_HAVE_GSL
  std::cout << "Eucken (with KineticsTheoryViscosity):" << std::endl;
  for( unsigned int s = 0; s < n_species; s++ )
    {
      Scalar mu = k_mu_mixture(s,T);
      Scalar value = e_k_mixture.conductivity_without_diffusion(s, T, mu);
      if( s == N2_species )
        k_e_kt_N2_value = value;

      std::cout << "k(" << species_str_list[s] << ") = " << value << std::endl;
    }
#endif // ANTIOCH_HAVE_GSL

  int return_flag = 0;

  // Now test all the values
  Scalar tol = std::numeric_limits<Scalar>::epsilon();

  if( !test_val<Scalar>(k_e_s_N2_value, k_e_s_N2_gold, tol, "Eucken-Sutherland") )
    return_flag = 1;

  if( !test_val<Scalar>(k_e_b_N2_value, k_e_b_N2_gold, tol, "Eucken-Blottner") )
    return_flag = 1;

#ifdef  ANTIOCH_HAVE_GSL
  if( !test_val<Scalar>(k_e_kt_N2_value, k_e_kt_N2_gold, tol, "Eucken-KineticsTheory") )
    return_flag = 1;
#endif // ANTIOCH_HAVE_GSL

  return return_flag;
}

int main()
{
  return (tester<double>() ||
          tester<long double>() ||
          tester<float>());
}
