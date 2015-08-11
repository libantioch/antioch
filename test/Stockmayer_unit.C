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
int check_value(const Scalar ref, const Scalar candidate, const Scalar tol, const std::string& test_name)
{
  int return_flag = 0;

  std::cout << "Testing "+test_name << "...";

  Scalar error = std::abs((ref - candidate)/ref);

  if( error > tol)
    {
      std::cout << "FAILED!" << std::endl;
      std::cerr << std::scientific << std::setprecision(16);
      std::cerr << "  reference = " << ref << std::endl
                << "  candidate = " << candidate << std::endl
                << "  relative error = " << error << std::endl
                << "  tolerance = " << tol << std::endl;
      return_flag = 1;
    }
  else
    {
      std::cout << "PASSED!" << std::endl;
    }

  return return_flag;
}

template <typename Scalar>
int tester()
{
  std::vector<std::string> molecules;
  molecules.push_back("CH4"); // T_max = 14,140 K; 11,743 K; 7,330 K
  molecules.push_back("N2");  // T_max = 9,753 K;  11,743 K; 6,088 K
  molecules.push_back("H2");  // T_max = 3,800 K;  6,8088 K; 7,330 K

  const unsigned int n_species = molecules.size();

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

  Scalar T_max(8000.L);

  ps_mu.extrapolate_max_temp(T_max);
  bimol_D.extrapolate_max_temp(T_max);

  std::vector<std::vector<Scalar> > D_matrix_reg(n_species);
  std::vector<std::vector<Scalar> > D_matrix(n_species);
  for( unsigned int s = 0; s < n_species; s++ )
    {
      D_matrix[s].resize(n_species);
      D_matrix_reg[s].resize(n_species);
    }

  Scalar mu_CH4_reg = 1.0573148339577483e-04;
  Scalar mu_N2_reg = 1.5770745584335467e-04;
  Scalar mu_H2_reg = 4.6078198688681625e-05;

  D_matrix_reg[0][0] = 8.1850066826705467e-03;
  D_matrix_reg[0][1] = 7.7255676543780075e-03;
  D_matrix_reg[0][2] = 4.5830263450987500e-02;
  D_matrix_reg[1][0] = 7.7255676543780075e-03;
  D_matrix_reg[1][1] = 7.0195276250601020e-03;
  D_matrix_reg[1][2] = 3.1562763815500983e-03;
  D_matrix_reg[2][0] = 4.5830263450987500e-02;
  D_matrix_reg[2][1] = 3.1562763815500983e-03;
  D_matrix_reg[2][2] = -4.3151204734397887e-03;

  Scalar mu_CH4 = ps_mu(0,7900.0);
  Scalar mu_N2 = ps_mu(1,7900.0);
  Scalar mu_H2 = ps_mu(2,3900.0);

  bimol_D.compute_binary_diffusion_matrix(6900.0, 1.0, D_matrix);

  Scalar tol = 10.0*std::numeric_limits<Scalar>::epsilon();

  // Check viscosity against regression values
  return_flag = check_value( mu_CH4_reg, mu_CH4, tol, "mu_CH4" ) ||
                check_value( mu_N2_reg, mu_N2, tol, "mu_N2" ) ||
                check_value( mu_H2_reg, mu_H2, tol, "mu_H2" );

  // Matrix had better be symmetric
  for( unsigned int s = 0; s < n_species; s++ )
    {
      for( unsigned int t = s+1; t < n_species; t++ )
        {
          std::stringstream convert;
          convert << s;
          std::string test_name = "D["+convert.str()+"]";
          convert.str("");
          convert << t;
          test_name += "["+convert.str()+"] symmetry";

          return_flag = return_flag ||
            check_value( D_matrix[s][t], D_matrix[t][s], tol, test_name );
        }
    }

  // Check binary diffusion matrix against regression values
  for( unsigned int s = 0; s < n_species; s++ )
    {
      for( unsigned int t = 0; t < n_species; t++ )
        {
          std::stringstream convert;
          convert << s;
          std::string test_name = "D["+convert.str()+"]";
          convert.str("");
          convert << t;
          test_name += "["+convert.str()+"]";

          return_flag = return_flag ||
            check_value( D_matrix_reg[s][t], D_matrix[s][t], tol, test_name );
        }
    }

  return return_flag;
}



#endif // ANTIOCH_HAVE_GSL


int main()
{
#ifdef ANTIOCH_HAVE_GSL
  // No long double test because GSL (used for splining) only supports double
  return (tester<float>() ||
          tester<double>());
#else
  // 77 return code tells Automake we skipped this.
  return 77;
#endif
}
