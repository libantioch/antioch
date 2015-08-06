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


#include "antioch_config.h"

// Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/mixture_viscosity.h"
#include "antioch/mixture_diffusion.h"
#include "antioch/molecular_binary_diffusion.h"
#include "antioch/kinetics_theory_viscosity.h"
#include "antioch/kinetics_theory_viscosity_building.h"
#include "antioch/cea_curve_fit.h"
#include "antioch/nasa_mixture.h"
#include "antioch/nasa_evaluator.h"
#include "antioch/ideal_gas_micro_thermo.h"
#include "antioch/transport_mixture.h"
#include "antioch/gsl_spliner.h"
#include "antioch/cea_mixture_parsing.h"
#include "antioch/default_filename.h"
#include "antioch/vector_utils.h"

// C++
#include <limits>
#include <iostream>
#include <iomanip>


#ifdef ANTIOCH_HAVE_GSL

// not necessary (yet?)
template <typename Scalar>
int check_value(const Scalar & ref, const Scalar & candidate, const Scalar & x, const std::string & words)
{
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 10.;

  if(std::abs((ref - candidate)/ref) > tol)
  {
      std::cerr << std::scientific << std::setprecision(15);
      std::cerr << "  reference = " << ref << std::endl
                << "  candidate = " << candidate << std::endl
                << "  relative difference = " << std::abs((ref - candidate) / ref) << std::endl
                << "  absolute difference = " << std::abs(ref - candidate) << std::endl
                << "  tolerance = " << tol << std::endl;
      return 1;
  }

  return 0;
}

template <typename Scalar>
int tester()
{
  std::vector<std::string> molecules;
  molecules.push_back("CH4"); // T_max = 14,140 K; 11,743 K; 7,330 K
  molecules.push_back("N2");  // T_max = 9,753 K;  11,743 K; 6,088 K
  molecules.push_back("H2");  // T_max = 3,800 K;  6,8088 K; 7,330 K

  // mixture
  Antioch::ChemicalMixture<Scalar> chem_mixture( molecules );

  // macro thermo
  Antioch::NASAThermoMixture<Scalar, Antioch::NASA9CurveFit<Scalar> > nasa_mixture( chem_mixture );

  // ASCII and true are default, but let's be verbose
  Antioch::read_nasa_mixture_data( nasa_mixture, Antioch::DefaultFilename::thermo_data(), Antioch::ASCII, true);

  typedef Antioch::NASAEvaluator<Scalar, Antioch::NASA9CurveFit<Scalar> > MacroThermo;
  MacroThermo nasa_thermo( nasa_mixture );

  // micro thermo
  typedef Antioch::IdealGasMicroThermo<MacroThermo,Scalar> MicroThermo;
  MicroThermo micro_thermo(nasa_thermo, chem_mixture);

  // transport
  Antioch::TransportMixture<Scalar> tran_mixture( chem_mixture );


  // sets, GSLSpliner implicit
  Antioch::MixtureViscosity<Antioch::KineticsTheoryViscosity<Scalar,Antioch::GSLSpliner>,Scalar >
    ps_mu(tran_mixture);
  Antioch::build_kinetics_theory_viscosity<Scalar,Antioch::GSLSpliner>(ps_mu);

  Antioch::MixtureDiffusion<Antioch::MolecularBinaryDiffusion<Scalar,Antioch::GSLSpliner>,Scalar>
    bimol_D( tran_mixture );

  int return_flag(0);

  //Scalar T_max(8000.L);

  //extrapolate_to_high_temperatures(viscosity,T_max);
  //extrapolate_to_high_temperatures(diffusion,T_max);

  return return_flag;
}



#endif // ANTIOCH_HAVE_GSL


int main()
{
#ifdef ANTIOCH_HAVE_GSL
// gsl work in double...
  return (tester<float>() ||
          tester<double>());
 //         tester<long double>() ||
#else
  // 77 return code tells Automake we skipped this.
  return 77;
#endif
}
