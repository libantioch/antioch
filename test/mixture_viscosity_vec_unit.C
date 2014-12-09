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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

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
#include "antioch/vexcl_utils_decl.h"

#include "antioch/default_filename.h"
#include "antioch/sutherland_viscosity.h"
#include "antioch/blottner_viscosity.h"
#include "antioch/blottner_viscosity_utils.h"
#include "antioch/sutherland_viscosity_utils.h"
#include "antioch/physics_metaprogramming.h"
#include "antioch/physical_set.h"
#include "antioch/blottner_parsing.h"
#include "antioch/sutherland_parsing.h"

#include "antioch/eigen_utils.h"
#include "antioch/metaphysicl_utils.h"
#include "antioch/valarray_utils.h"
#include "antioch/vexcl_utils.h"

#ifdef ANTIOCH_HAVE_GRVY
#include "grvy.h"

GRVY::GRVY_Timer_Class gt;
#endif

// C++
#include <cmath>
#include <iostream>

template <typename PairScalars>
int vectester(const PairScalars& example, const std::string& testname)
{
  typedef typename Antioch::value_type<PairScalars>::type Scalar;

  std::vector<std::string> species_str_list;
  const unsigned int n_species = 2;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "Air" ); // Yes, I know this doesn't make sense, it's just a test.

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );

  Antioch::PhysicalSet<Antioch::SutherlandViscosity<Scalar>, Antioch::ChemicalMixture<Scalar> > s_mu_mixture(chem_mixture);

  Antioch::PhysicalSet<Antioch::BlottnerViscosity<Scalar>, Antioch::ChemicalMixture<Scalar> > b_mu_mixture(chem_mixture);

  Antioch::read_sutherland_data_ascii<Scalar>( s_mu_mixture, Antioch::DefaultFilename::sutherland_data() );
  Antioch::read_blottner_data_ascii<Scalar>( b_mu_mixture, Antioch::DefaultFilename::blottner_data() );

  std::cout << s_mu_mixture << std::endl;
  std::cout << b_mu_mixture << std::endl;

  PairScalars T = example;
  for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
    {
      T[2*tuple  ] = 1500.1;
      T[2*tuple+1] = 1600.1;
    }

  PairScalars mu = example;

  std::cout << "Blottner:" << std::endl;
  for( unsigned int s = 0; s < n_species; s++ )
    {
#ifdef ANTIOCH_HAVE_GRVY
      const std::string testblottner = testname + "-blottner";
      gt.BeginTimer(testblottner);
#endif

      b_mu_mixture(s,T,mu);

#ifdef ANTIOCH_HAVE_GRVY
      gt.EndTimer(testblottner);
#endif

      std::cout << "mu(" << species_str_list[s] << ") = " << mu << std::endl;
    }

  std::cout << "Sutherland:" << std::endl;
  for( unsigned int s = 0; s < n_species; s++ )
    {
#ifdef ANTIOCH_HAVE_GRVY
      const std::string testsutherland = testname + "-sutherland";
      gt.BeginTimer(testsutherland);
#endif

      s_mu_mixture(s,T,mu);

#ifdef ANTIOCH_HAVE_GRVY
      gt.EndTimer(testsutherland);
#endif

      std::cout << "mu(" << species_str_list[s] << ") = " << mu << std::endl;
    }

  int return_flag = 0;

  return return_flag;
}

int main()
{
  int returnval = 0;

  returnval = returnval ||
    vectester (std::valarray<float>(2*ANTIOCH_N_TUPLES), "valarray<float>");
  returnval = returnval ||
    vectester (std::valarray<double>(2*ANTIOCH_N_TUPLES), "valarray<double>");
  returnval = returnval ||
    vectester (std::valarray<long double>(2*ANTIOCH_N_TUPLES), "valarray<ld>");
#ifdef ANTIOCH_HAVE_EIGEN
  returnval = returnval ||
    vectester (Eigen::Array<float, 2*ANTIOCH_N_TUPLES, 1>(), "Eigen::ArrayXf");
  returnval = returnval ||
    vectester (Eigen::Array<double, 2*ANTIOCH_N_TUPLES, 1>(), "Eigen::ArrayXd");
  returnval = returnval ||
    vectester (Eigen::Array<long double, 2*ANTIOCH_N_TUPLES, 1>(), "Eigen::ArrayXld");
#endif
#ifdef ANTIOCH_HAVE_METAPHYSICL
  returnval = returnval ||
    vectester (MetaPhysicL::NumberArray<2*ANTIOCH_N_TUPLES, float> (0), "NumberArray<float>");
  returnval = returnval ||
    vectester (MetaPhysicL::NumberArray<2*ANTIOCH_N_TUPLES, double> (0), "NumberArray<double>");
  returnval = returnval ||
    vectester (MetaPhysicL::NumberArray<2*ANTIOCH_N_TUPLES, long double> (0), "NumberArray<ld>");
#endif
#ifdef ANTIOCH_HAVE_VEXCL
  vex::Context ctx_f (vex::Filter::All);
  if (!ctx_f.empty())
    returnval = returnval ||
      vectester (vex::vector<float> (ctx_f, 2*ANTIOCH_N_TUPLES), "vex::vector<float>");

  vex::Context ctx_d (vex::Filter::DoublePrecision);
  if (!ctx_d.empty())
    returnval = returnval ||
      vectester (vex::vector<double> (ctx_d, 2*ANTIOCH_N_TUPLES), "vex::vector<double>");
#endif

#ifdef ANTIOCH_HAVE_GRVY
  gt.Finalize();
  gt.Summarize();
#endif

  return returnval;
}
