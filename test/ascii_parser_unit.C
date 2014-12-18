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

#include "antioch/chemical_mixture.h"
#include "antioch/transport_mixture.h"
#include "antioch/ascii_parser.h"
#include "antioch/ideal_gas_internal_thermo.h"
#include "antioch/cea_curve_fit.h"
#include "antioch/nasa_mixture.h"
#include "antioch/nasa_evaluator.h"
#include "antioch/thermo_evaluator.h"

#include <iomanip>
#include <limits>
#include <string>

template <typename Scalar>
int test_value(const Scalar & value, const Scalar & ref, const std::string & words)
{
   const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 10;
   const Scalar rel_diff = (ref < tol)?value:
                                       std::abs((ref - value)/ref);

   if(rel_diff > tol)
   {
      std::cerr << std::setprecision(16)
                << "Error in comparison of " << words    << "\n"
                << "Reference value = "      << ref      << "\n"
                << "Computed value = "       << value    << "\n"
                << "Relative distance = "    << rel_diff << "\n"
                << "Tolerance = "            << tol      << std::endl;
      return 1;
   }


   return 0;
}

template <typename Scalar>
int test_transport_species(const Antioch::TransportSpecies<Scalar> & spec, Scalar depth, Scalar sigma, Scalar dip, Scalar pol, Scalar Zrot)
{
   int return_flag(0);

   return_flag = test_value(spec.LJ_depth(),depth, "Lennard-Jones depth")               || return_flag;
   return_flag = test_value(spec.LJ_diameter(),sigma, "Lennard-Jones radius")           || return_flag;
   return_flag = test_value(spec.dipole_moment(),dip, "Dipole moment")                  || return_flag;
   return_flag = test_value(spec.polarizability(),pol, "polarizability")                || return_flag;
   return_flag = test_value(spec.rotational_relaxation(),Zrot, "Rotational relaxation") || return_flag;

   return return_flag;
}

template <typename Scalar>
int tester(const std::string& filename)
{

   std::vector<std::string> species_list;
   species_list.push_back("N2");
   species_list.push_back("O2");
   species_list.push_back("H2");

    Antioch::ChemicalMixture<Scalar> chem_mix(species_list);
//thermo
//// macro
    Antioch::NASAThermoMixture<Scalar, Antioch::CEACurveFit<Scalar> > nasa_mixture( chem_mix);
    Antioch::NASAEvaluator<Scalar, Antioch::CEACurveFit<Scalar> >     nasa_thermo( nasa_mixture );

    typedef Antioch::NASAEvaluator<Scalar, Antioch::CEACurveFit<Scalar> > NASAThermoType;

//// micro
    Antioch::IdealGasInternalThermo<NASAThermoType,Scalar>  micro_thermo(nasa_thermo,chem_mix);

    typedef Antioch::IdealGasInternalThermo<NASAThermoType,Scalar>  StatThermoType;

//// eval
    Antioch::ThermoEvaluator<Scalar,NASAThermoType,StatThermoType> thermo(nasa_thermo,micro_thermo);

    typedef Antioch::ThermoEvaluator<Scalar,NASAThermoType,StatThermoType> ThermoEval;

//transport
///// mixture and stat thermo for conduction
        // read
    Antioch::ASCIIParser<Scalar> parser(filename,true);

    parser.set_ignored_columns(std::vector<unsigned int>(1,0));

    Antioch::TransportMixture<ThermoEval,Scalar> tran_mixture( chem_mix, thermo, parser);
  
   int return_flag(0);
// name          geom      LJ_depth  LJ_sigma   dipole   polar          Zrot
    // testing N2
// N2                 1    97.530     3.621     0.000     1.760     4.000

    return_flag = test_transport_species(*(tran_mixture.transport_species()[0]), 97.530, 3.621, 0.000, 1.760, 4.000) || return_flag;

    // testing O2
// O2                 1   107.400     3.458     0.000     1.600     3.800
    return_flag = test_transport_species(*(tran_mixture.transport_species()[1]), 107.400, 3.458, 0.000, 1.600, 3.800) || return_flag;

    // testing H2
// H2                 1    38.000     2.920     0.000     0.790   280.000
    return_flag = test_transport_species(*(tran_mixture.transport_species()[2]), 38.000, 2.920, 0.000, 0.790, 280.000) || return_flag;

    return return_flag;
}

int main(int argc, char* argv[])
{
  // Check command line count.
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify ascii transport input file." << std::endl;
      antioch_error();
    }

  return tester<double>(std::string(argv[1]));
}
