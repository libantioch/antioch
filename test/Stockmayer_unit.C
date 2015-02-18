//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-


#include "antioch_config.h"

// Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/physics_metaprogramming_decl.h"
#include "antioch/physics_utils.h"
#include "antioch/molecular_binary_diffusion.h"
#include "antioch/kinetics_theory_viscosity.h"
#include "antioch/cea_curve_fit.h"
#include "antioch/cea_mixture.h"
#include "antioch/cea_evaluator.h"
#include "antioch/ideal_gas_micro_thermo.h"
#include "antioch/thermo_handler.h"
#include "antioch/transport_mixture.h"
#include "antioch/gsl_spliner.h"
#include "antioch/physical_set.h"
#include "antioch/cea_mixture_ascii_parsing.h"
#include "antioch/kinetics_theory_viscosity_building.h"
#include "antioch/molecular_binary_diffusion_building.h"
#include "antioch/default_filename.h"
#include "antioch/physics_metaprogramming.h"
#include "antioch/kinetics_theory_viscosity_utils.h"
#include "antioch/molecular_binary_diffusion_utils.h"
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
  Antioch::CEAThermoMixture<Scalar> nasa_mixture( chem_mixture );
  Antioch::read_cea_mixture_data_ascii( nasa_mixture, Antioch::DefaultFilename::thermo_data());
  Antioch::CEAEvaluator<Scalar > nasa_thermo( nasa_mixture );

  typedef Antioch::CEAEvaluator<Scalar> MacroThermo;

  // micro thermo 
  Antioch::IdealGasMicroThermo<MacroThermo,Scalar> micro_thermo(nasa_thermo, chem_mixture);

  typedef Antioch::IdealGasMicroThermo<MacroThermo,Scalar> MicroThermo;

  // thermo, micro is implicit 
  Antioch::ThermoHandler<Scalar,MacroThermo > thermo(nasa_thermo,micro_thermo);

  typedef Antioch::ThermoHandler<Scalar,MacroThermo, MicroThermo > Thermo;

  // transport
  Antioch::TransportMixture<Thermo,Scalar> transport(chem_mixture, thermo);
  Antioch::read_transport_species_data_ascii(transport, Antioch::DefaultFilename::transport_mixture());

  typedef Antioch::TransportMixture<Thermo,Scalar> Transport;

  // sets, GSLSpliner implicit
  Antioch::PhysicalSet<Antioch::KineticsTheoryViscosity<Scalar>,Transport> viscosity(transport);

  Antioch::PhysicalSet<Antioch::MolecularBinaryDiffusion<Scalar>,Transport> diffusion(transport);

  int return_flag(0);

  Scalar T_max(8000.L);
  
  extrapolate_to_high_temperatures(viscosity,T_max);
  extrapolate_to_high_temperatures(diffusion,T_max);

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
