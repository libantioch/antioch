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

#include "antioch_config.h"

#include <valarray>

#ifdef ANTIOCH_HAVE_EIGEN
#include "Eigen/Dense"
#endif

#ifdef ANTIOCH_HAVE_METAPHYSICL
#include "metaphysicl/numberarray.h"
#endif

#ifdef ANTIOCH_HAVE_VEXCL
#include "vexcl/vexcl.hpp"
#endif

// Antioch

// Declare metaprogramming overloads before they're used
#include "antioch/eigen_utils_decl.h"
#include "antioch/metaphysicl_utils_decl.h"
#include "antioch/valarray_utils_decl.h"
#include "antioch/vector_utils_decl.h"
#include "antioch/vexcl_utils_decl.h"

#include "antioch/antioch_asserts.h"
#include "antioch/chemical_species.h"
#include "antioch/chemical_mixture.h"
#include "antioch/reaction_set.h"
#include "antioch/read_reaction_set_data.h"
#include "antioch/cea_evaluator.h"
#include "antioch/cea_mixture_ascii_parsing.h"
#include "antioch/kinetics_evaluator.h"

#include "antioch/eigen_utils.h"
#include "antioch/metaphysicl_utils.h"
#include "antioch/valarray_utils.h"
#include "antioch/vector_utils.h"
#include "antioch/vexcl_utils.h"

#ifdef ANTIOCH_HAVE_GRVY
#include "grvy.h"

GRVY::GRVY_Timer_Class gt;
#endif

// C++
#include <limits>
#include <string>
#include <vector>

static const unsigned int n_species = 5;

template <typename SpeciesVector1, typename SpeciesVector2>
int vec_compare (const SpeciesVector1 &a, const SpeciesVector2 &b, const std::string &name)
{
  int found_errors = 0;

  typedef typename Antioch::value_type<SpeciesVector1>::type StateType;
  typedef typename Antioch::value_type<StateType>::type Scalar;

  if (static_cast<std::size_t>(a.size()) !=
      static_cast<std::size_t>(b.size()))
    {
      std::cerr << "Error: Mismatch in vector sizes " << name << std::endl;
    }

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 200;

  for (unsigned int s=0; s != a.size(); s++)
    {
      using std::abs;
      using std::max;

      // Break this expression up to workaround bugs in my Eigen
      // and VexCL versions - RHS
      const StateType as = a[s], bs = b[s];
      const StateType rel_error = (as - bs)/max(as,bs);
      const StateType abs_rel_error = abs(rel_error);

      if( Antioch::max(abs_rel_error) > tol )
	{
	  found_errors++;
	}
    }

  if (found_errors)
    {
      using std::max;

      std::cerr << "Error: Mismatch in vectors " << name << std::endl;
      for( unsigned int s = 0; s < n_species; s++)
	{
          const StateType as = a[s], bs = b[s];
          const StateType rel_error = (as - bs)/max(as,bs);
	  std::cout << std::scientific << std::setprecision(16)
		    << "a(" << s << ") = " << as << std::endl
		    << "b(" << s << ") = " << bs << std::endl
		    << "a-b(" << s << ") = " << StateType(as-bs)
		    << std::endl <<
		    "a-b/(max(a,b)) = " << rel_error << std::endl
		    << "tol =" << tol
		    << std::endl;
	}
    }

  return found_errors;
}

template <typename PairScalars>
int vectester(const std::string& input_name,
	      const PairScalars& example,
	      const std::string& testname )
{
  typedef typename Antioch::value_type<PairScalars>::type Scalar;

  std::vector<std::string> species_str_list;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "O2" );
  species_str_list.push_back( "N" );
  species_str_list.push_back( "O" );
  species_str_list.push_back( "NO" );

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );
  Antioch::ReactionSet<Scalar> reaction_set( chem_mixture );

  Antioch::CEAThermoMixture<Scalar> cea_mixture( chem_mixture );
  Antioch::read_cea_mixture_data_ascii( cea_mixture, Antioch::DefaultFilename::thermo_data() );
  Antioch::CEAEvaluator<Scalar> thermo( cea_mixture );

  Antioch::read_reaction_set_data_xml<Scalar>( input_name, true, reaction_set );

  PairScalars T = example;
  PairScalars P = example;

  // Mass fractions
  PairScalars massfrac = example;

  for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
    {
      T[2*tuple  ] = 1500.0;
      T[2*tuple+1] = 1500.0;

      P[2*tuple  ] = 1.0e5;
      P[2*tuple+1] = 1.0e5;

      massfrac[2*tuple  ] = 0.2;
      massfrac[2*tuple+1] = 0.2;
    }

  const Antioch::KineticsConditions<PairScalars> conditions(T);
  const std::vector<PairScalars> Y(n_species,massfrac);
  std::vector<PairScalars> molar_densities(n_species, example);
  std::vector<PairScalars> h_RT_minus_s_R(n_species, example);
  std::vector<PairScalars> dh_RT_minus_s_R_dT(n_species, example);

  Antioch::KineticsEvaluator<Scalar,PairScalars> kinetics( reaction_set, example );

  std::vector<PairScalars> omega_dot(n_species, example);
  std::vector<PairScalars> omega_dot_2(n_species, example);
  std::vector<PairScalars> domega_dot_dT(n_species, example);

  std::vector<std::vector<PairScalars> > domega_dot_drhos
    (n_species, omega_dot); // omega_dot is a good example vec<Pair>


#ifdef ANTIOCH_HAVE_GRVY
  const std::string testnormal = testname + "-normal";
  gt.BeginTimer(testnormal);
#endif

  const PairScalars R_mix = chem_mixture.R(Y);

  const PairScalars rho = P/(R_mix*T);

  chem_mixture.molar_densities(rho,Y,molar_densities);

  Antioch::TempCache<PairScalars> temp_cache(T);
  thermo.h_RT_minus_s_R(temp_cache,h_RT_minus_s_R);
  thermo.dh_RT_minus_s_R_dT(temp_cache,dh_RT_minus_s_R_dT);

  kinetics.compute_mass_sources( conditions, molar_densities, h_RT_minus_s_R, omega_dot );

  kinetics.compute_mass_sources_and_derivs ( conditions,
					     molar_densities,
					     h_RT_minus_s_R,
					     dh_RT_minus_s_R_dT,
					     omega_dot_2,
					     domega_dot_dT,
					     domega_dot_drhos );

#ifdef ANTIOCH_HAVE_GRVY
  gt.EndTimer(testnormal);
#endif

  int return_flag = 0;

#ifdef ANTIOCH_HAVE_EIGEN
  {
    typedef Eigen::Array<PairScalars,n_species,1> SpeciesVecEigenType;

#ifdef ANTIOCH_HAVE_GRVY
    const std::string testeigenA = testname + "-eigenA";
    gt.BeginTimer(testeigenA);
#endif

    const Antioch::KineticsConditions<PairScalars,SpeciesVecEigenType> eigen_conditions(T);

    SpeciesVecEigenType eigen_Y;
    Antioch::init_constant(eigen_Y, massfrac);

    SpeciesVecEigenType eigen_molar_densities;
    Antioch::init_constant(eigen_molar_densities, example);

    SpeciesVecEigenType eigen_h_RT_minus_s_R;
    Antioch::init_constant(eigen_h_RT_minus_s_R, example);

    SpeciesVecEigenType eigen_dh_RT_minus_s_R_dT;
    Antioch::init_constant(eigen_dh_RT_minus_s_R_dT, example);

    SpeciesVecEigenType eigen_omega_dot;
    Antioch::init_constant(eigen_omega_dot, example);

    SpeciesVecEigenType eigen_domega_dot_dT;
    Antioch::init_constant(eigen_domega_dot_dT, example);

    // FIXME: What to do for domega_dot_drhos type?

    const PairScalars eigen_R = chem_mixture.R(eigen_Y);
    chem_mixture.molar_densities(rho,eigen_Y,eigen_molar_densities);

    thermo.h_RT_minus_s_R(temp_cache,eigen_h_RT_minus_s_R);

    thermo.dh_RT_minus_s_R_dT(temp_cache,eigen_dh_RT_minus_s_R_dT);

    kinetics.compute_mass_sources( eigen_conditions, eigen_molar_densities, eigen_h_RT_minus_s_R, eigen_omega_dot );

#ifdef ANTIOCH_HAVE_GRVY
    gt.EndTimer(testeigenA);
#endif

    return_flag +=
      vec_compare(eigen_R,R_mix,"eigen_R");

    return_flag +=
      vec_compare(eigen_molar_densities, molar_densities,
                  "eigen_molar_densities");

    return_flag +=
      vec_compare(eigen_h_RT_minus_s_R, h_RT_minus_s_R,
                  "eigen_h_RT_minus_s_R");

    return_flag +=
      vec_compare(eigen_dh_RT_minus_s_R_dT, dh_RT_minus_s_R_dT,
                  "eigen_dh_RT_minus_s_R_dT");

    return_flag +=
      vec_compare(eigen_omega_dot,omega_dot,"eigen_omega_dot");
  }

  {
    typedef Eigen::Matrix<PairScalars,Eigen::Dynamic,1> SpeciesVecEigenType;
    const Antioch::KineticsConditions<PairScalars,SpeciesVecEigenType> eigen_conditions(T);

    SpeciesVecEigenType eigen_Y(n_species,1);
    Antioch::init_constant(eigen_Y, massfrac);

    SpeciesVecEigenType eigen_molar_densities(n_species,1);
    Antioch::init_constant(eigen_molar_densities, example);

    SpeciesVecEigenType eigen_h_RT_minus_s_R(n_species,1);
    Antioch::init_constant(eigen_h_RT_minus_s_R, example);

    SpeciesVecEigenType eigen_dh_RT_minus_s_R_dT(n_species,1);
    Antioch::init_constant(eigen_dh_RT_minus_s_R_dT, example);

    SpeciesVecEigenType eigen_omega_dot(n_species,1);
    Antioch::init_constant(eigen_omega_dot, example);

    SpeciesVecEigenType eigen_domega_dot_dT(n_species,1);
    Antioch::init_constant(eigen_domega_dot_dT, example);

    // FIXME: What to do for domega_dot_drhos type?

#ifdef ANTIOCH_HAVE_GRVY
    const std::string testeigenV = testname + "-eigenV";
    gt.BeginTimer(testeigenV);
#endif

    const PairScalars eigen_R = chem_mixture.R(eigen_Y);

    chem_mixture.molar_densities(rho,eigen_Y,eigen_molar_densities);

    thermo.h_RT_minus_s_R(temp_cache,eigen_h_RT_minus_s_R);

    thermo.dh_RT_minus_s_R_dT(temp_cache,eigen_dh_RT_minus_s_R_dT);

    kinetics.compute_mass_sources( eigen_conditions, eigen_molar_densities, eigen_h_RT_minus_s_R, eigen_omega_dot );

#ifdef ANTIOCH_HAVE_GRVY
    gt.EndTimer(testeigenV);
#endif

    return_flag +=
      vec_compare(eigen_R,R_mix,"eigen_R");

    return_flag +=
      vec_compare(eigen_molar_densities, molar_densities,
                  "eigen_molar_densities");

    return_flag +=
      vec_compare(eigen_h_RT_minus_s_R, h_RT_minus_s_R,
                  "eigen_h_RT_minus_s_R");

    return_flag +=
      vec_compare(eigen_dh_RT_minus_s_R_dT, dh_RT_minus_s_R_dT,
                  "eigen_dh_RT_minus_s_R_dT");

    return_flag +=
      vec_compare(eigen_omega_dot,omega_dot,"eigen_omega_dot");
  }

#endif // ANTIOCH_HAVE_EIGEN

  for( unsigned int s = 0; s < n_species; s++)
    {
      std::cout << std::scientific << std::setprecision(16)
		<< "omega_dot(" << chem_mixture.chemical_species()[s]->species() << ") = "
		<< omega_dot[s] << std::endl;
    }

  for( unsigned int s = 0; s < n_species; s++)
    {
      std::cout << std::scientific << std::setprecision(16)
		<< "domega_dot_dT(" << chem_mixture.chemical_species()[s]->species() << ") = "
		<< domega_dot_dT[s] << std::endl;
    }

  for( unsigned int s = 0; s < n_species; s++)
    {
      for( unsigned int t = 0; t < n_species; t++)
        {
          std::cout << std::scientific << std::setprecision(16)
		    << "domega_dot_drhos(" << chem_mixture.chemical_species()[s]->species()
		    << ", " << chem_mixture.chemical_species()[t]->species() << ") = "
		    << domega_dot_drhos[s][t] << std::endl;
        }
    }

  // Regression values for omega_dot
  std::vector<PairScalars> omega_dot_reg(n_species,example);

  // Regression values for domega_dot_drhos
  std::vector<std::vector<PairScalars> > domega_dot_drhos_reg
    (n_species, omega_dot); // omega_dot is the right size for an example

  // Regression values for domega_dot_dT
  std::vector<PairScalars> domega_dot_dT_reg(n_species, example);

  for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
    {
      omega_dot_reg[0][2*tuple  ] =  9.1623705357123753e+04;
      omega_dot_reg[1][2*tuple  ] = -3.3462025680272243e+05;
      omega_dot_reg[2][2*tuple  ] = -2.1139216712069495e+05;
      omega_dot_reg[3][2*tuple  ] =  1.9782018625609628e+05;
      omega_dot_reg[4][2*tuple  ] =  2.5656853231019735e+05;
      omega_dot_reg[0][2*tuple+1] = Scalar(omega_dot_reg[0][0]);
      omega_dot_reg[1][2*tuple+1] = Scalar(omega_dot_reg[1][0]);
      omega_dot_reg[2][2*tuple+1] = Scalar(omega_dot_reg[2][0]);
      omega_dot_reg[3][2*tuple+1] = Scalar(omega_dot_reg[3][0]);
      omega_dot_reg[4][2*tuple+1] = Scalar(omega_dot_reg[4][0]);

      domega_dot_dT_reg[0][2*tuple  ] =  1.8014990183270937e+02;
      domega_dot_dT_reg[1][2*tuple  ] = -5.2724437115534380e+02;
      domega_dot_dT_reg[2][2*tuple  ] = -3.0930094476883017e+02;
      domega_dot_dT_reg[3][2*tuple  ] =  3.7972747459781005e+02;
      domega_dot_dT_reg[4][2*tuple  ] =  2.7666793949365456e+02;
      domega_dot_dT_reg[0][2*tuple+1] = Scalar(domega_dot_dT_reg[0][0]);
      domega_dot_dT_reg[1][2*tuple+1] = Scalar(domega_dot_dT_reg[1][0]);
      domega_dot_dT_reg[2][2*tuple+1] = Scalar(domega_dot_dT_reg[2][0]);
      domega_dot_dT_reg[3][2*tuple+1] = Scalar(domega_dot_dT_reg[3][0]);
      domega_dot_dT_reg[4][2*tuple+1] = Scalar(domega_dot_dT_reg[4][0]);

      domega_dot_drhos_reg[0][0][2*tuple  ] = 1.9675775188085109e+04;
      domega_dot_drhos_reg[0][1][2*tuple  ] = 1.7226141262419737e+04;
      domega_dot_drhos_reg[0][2][2*tuple  ] = 3.2159299284723610e+06;
      domega_dot_drhos_reg[0][3][2*tuple  ] = 1.4765214711933021e+05;
      domega_dot_drhos_reg[0][4][2*tuple  ] = 2.3225053279918131e+06;

      domega_dot_drhos_reg[1][0][2*tuple  ] =  8.8927385505978492e+03;
      domega_dot_drhos_reg[1][1][2*tuple  ] = -9.9560178070099482e+06;
      domega_dot_drhos_reg[1][2][2*tuple  ] = -9.8748760140991123e+06;
      domega_dot_drhos_reg[1][3][2*tuple  ] =  4.6143036700500813e+05;
      domega_dot_drhos_reg[1][4][2*tuple  ] =  8.3487375168772399e+03;

      domega_dot_drhos_reg[2][0][2*tuple  ] = -2.2420842426881281e+04;
      domega_dot_drhos_reg[2][1][2*tuple  ] = -4.3812843857644886e+06;
      domega_dot_drhos_reg[2][2][2*tuple  ] = -6.8343593463263955e+06;
      domega_dot_drhos_reg[2][3][2*tuple  ] = -5.4143671040862988e+05;
      domega_dot_drhos_reg[2][4][2*tuple  ] = -1.2267997668149246e+06;

      domega_dot_drhos_reg[3][0][2*tuple  ] = -1.2028166578920147e+04;
      domega_dot_drhos_reg[3][1][2*tuple  ] =  4.9713710400172938e+06;
      domega_dot_drhos_reg[3][2][2*tuple  ] =  5.7418898143800552e+06;
      domega_dot_drhos_reg[3][3][2*tuple  ] = -9.1121284934572734e+05;
      domega_dot_drhos_reg[3][4][2*tuple  ] =  1.2431710353864791e+06;

      domega_dot_drhos_reg[4][0][2*tuple  ] =  5.8804952671184686e+03;
      domega_dot_drhos_reg[4][1][2*tuple  ] =  9.3487050114947233e+06;
      domega_dot_drhos_reg[4][2][2*tuple  ] =  7.7514156175730915e+06;
      domega_dot_drhos_reg[4][3][2*tuple  ] =  8.4356704563001888e+05;
      domega_dot_drhos_reg[4][4][2*tuple  ] = -2.3472253340802449e+06;

      for (unsigned int i = 0; i != 5; ++i)
        for (unsigned int j = 0; j != 5; ++j)
          domega_dot_drhos_reg[i][j][2*tuple+1] =
            Scalar(domega_dot_drhos_reg[i][j][2*tuple  ]);
    }

  return_flag +=
    vec_compare(omega_dot, omega_dot_reg, "omega_dot");

  return_flag +=
    vec_compare(omega_dot_2, omega_dot_reg, "omega_dot_2");

  return_flag +=
    vec_compare(domega_dot_dT, domega_dot_dT_reg, "domega_dot_dT");

  char vecname[] = "domega_dot_drho_0";
  for (unsigned int i = 0; i != 5; ++i)
    {
      return_flag +=
        vec_compare(domega_dot_drhos[i], domega_dot_drhos_reg[i], vecname);

      // Should b=e safe to assume "1"+1 = "2" etc.
      ++vecname[16];
    }

  return return_flag;
}

int main(int argc, char* argv[])
{
  // Check command line count.
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify reaction set XML input file." << std::endl;
      antioch_error();
    }

  int returnval = 0;

  returnval +=
    vectester (argv[1], std::valarray<float>(2*ANTIOCH_N_TUPLES), "valarray<float>");
  returnval +=
    vectester (argv[1], std::valarray<double>(2*ANTIOCH_N_TUPLES), "valarray<double>");
// We're not getting the full long double precision yet?
//  returnval = returnval ||
//    vectester (argv[1], std::valarray<long double>(2*ANTIOCH_N_TUPLES), "valarray<ld>");
#ifdef ANTIOCH_HAVE_EIGEN
  returnval +=
    vectester (argv[1], Eigen::Array<float, 2*ANTIOCH_N_TUPLES, 1>(), "Eigen::ArrayXf");
  returnval +=
    vectester (argv[1], Eigen::Array<double, 2*ANTIOCH_N_TUPLES, 1>(), "Eigen::ArrayXd");
// We're not getting the full long double precision yet?
//  returnval = returnval ||
//    vectester (argv[1], Eigen::Array<long double, 2*ANTIOCH_N_TUPLES, 1>(), "Eigen::ArrayXld");
#endif
#ifdef ANTIOCH_HAVE_METAPHYSICL
  returnval +=
    vectester (argv[1], MetaPhysicL::NumberArray<2*ANTIOCH_N_TUPLES, float> (0), "NumberArray<float>");
  returnval +=
    vectester (argv[1], MetaPhysicL::NumberArray<2*ANTIOCH_N_TUPLES, double> (0), "NumberArray<double>");
//  returnval = returnval ||
//    vectester (argv[1], MetaPhysicL::NumberArray<2*ANTIOCH_N_TUPLES, long double> (0), "NumberArray<ld>");
#endif
#ifdef ANTIOCH_HAVE_VEXCL
  vex::Context ctx_f (vex::Filter::All);
  if (!ctx_f.empty())
    returnval = returnval ||
      vectester (argv[1], vex::vector<float> (ctx_f, 2*ANTIOCH_N_TUPLES), "vex::vector<float>");

  vex::Context ctx_d (vex::Filter::DoublePrecision);
  if (!ctx_d.empty())
    returnval = returnval ||
      vectester (argv[1], vex::vector<double> (ctx_d, 2*ANTIOCH_N_TUPLES), "vex::vector<double>");
#endif

  std::cout << "Found " << returnval << " errors" << std::endl;

#ifdef ANTIOCH_HAVE_GRVY
  gt.Finalize();
  gt.Summarize();
#endif

  return (returnval != 0) ? 1 : 0;
}
