//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

// C++
#include <cmath>
#include <limits>

// Antioch
#include "antioch/vector_utils_decl.h"

#include "antioch/physical_constants.h"
#include "antioch/chemical_mixture.h"
#include "antioch/ideal_gas_micro_thermo.h"
#include "antioch/cea_evaluator.h"
#include "antioch/cea_curve_fit.h"
#include "antioch/cea_mixture_ascii_parsing.h"

#include "antioch/vector_utils.h"

template <typename Scalar>
bool test_zero(const Scalar val, const Scalar tol)
{
  using std::abs;

  if( abs(val) > tol )
    return false;
  else
    return true;
}

template <typename Scalar>
int test_relative(const Scalar val, const Scalar truth, const Scalar tol, const std::string & words)
{
  using std::abs;

  bool test = (test_zero(truth,tol))?!test_zero(val,tol):(abs( (val-truth)/truth ) > tol );
  Scalar diff = (test_zero(truth,tol))?val:abs( (val-truth)/truth );
  if(test)
  {
    std::cerr << std::scientific << std::setprecision(20);
    std::cerr << "Error: Mismatch in " << words
              << "\n Expected      = " << truth
              << "\n Computed      = " << val
              << "\n Relative diff = " << diff
              << "\n Tolerance     = " << tol
              << std::endl;
    return 1;
  }else
  {
    return 0;
  }
}

template <typename Scalar>
int test_molecule(unsigned int spec, const Scalar & n_tr_dofs, const Scalar & R_spec, 
                   const Scalar & cv_over_R, const Scalar & T,
                   const Antioch::IdealGasMicroThermo< Antioch::NASAEvaluator<Scalar, Antioch::CEACurveFit<Scalar> >,
                                                          Scalar 
                                                        > & thermo,
                   const std::string & name)
{
    const Scalar tol = (std::numeric_limits<Scalar>::epsilon() * 2 < 5e-17L)?5e-17L:
                                                                             std::numeric_limits<Scalar>::epsilon() * 2;

    int return_flag = 0;
  // values over R are:
  // cv_tr / R    = n_tr_dofs
  // cv_trans / R = 1.5
  // cv_rot / R   = max(n_tr_dofs - 1.5,0)
  // cv_vib / R   = cp / R - 1 - n_tr_dofs

    Scalar cv_rot = (n_tr_dofs > Scalar(1.5))?n_tr_dofs - Scalar(1.5):Scalar(0.);
    Scalar cv_vib = (n_tr_dofs < 2.)?0.L:cv_over_R - n_tr_dofs;

    return_flag = test_relative(thermo.cv_tr_over_R(spec), n_tr_dofs, tol, "cv_tr_over_R of " + name)                  || return_flag;
    return_flag = test_relative(thermo.cv_tr(spec), R_spec * n_tr_dofs, tol, "cv_tr of " + name)                       || return_flag;
    return_flag = test_relative(thermo.cv_trans_over_R(spec), Scalar(1.5), tol, "cv_trans_over_R of " + name)          || return_flag;
    return_flag = test_relative(thermo.cv_trans(spec), R_spec * Scalar(1.5), tol, "cv_trans of " + name)               || return_flag;
    return_flag = test_relative(thermo.cv_rot_over_R(spec), cv_rot, tol, "cv_rot_over_R of " + name)                   || return_flag;
    return_flag = test_relative(thermo.cv_rot(spec), R_spec * cv_rot, tol, "cv_rot of " + name)                        || return_flag;
// vibration requires CEA fits, tolerance is somewhat loose...
    return_flag = test_relative(thermo.cv_vib_over_R(spec,T), cv_vib, Scalar(200.L) * tol, "cv_vib_over_R of " + name) || return_flag;
    return_flag = test_relative(thermo.cv_vib(spec,T),  R_spec * cv_vib, Scalar(200.L) * tol, "cv_vib of " + name)     || return_flag;

    return return_flag;
}
template <typename Scalar>
int tester()
{
  const Scalar Mm_N  = 14.008e-3;   //in SI kg/mol
  const Scalar Mm_O  = 16e-3;       //in SI kg/mol
  const Scalar Mm_N2 = 2.L * Mm_N;  //in SI kg/mol
  const Scalar Mm_O2 = 2.L * Mm_O;  //in SI kg/mol
  const Scalar Mm_NO = Mm_O + Mm_N; //in SI kg/mol

  std::vector<std::string> species_str_list;
  const unsigned int n_species = 5;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "O2" );
  species_str_list.push_back( "N" );
  species_str_list.push_back( "O" );
  species_str_list.push_back( "NO" );

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );

// required for Cv, we take default CEA
  Antioch::CEAThermoMixture<Scalar> nasa_mixture( chem_mixture );
  Antioch::read_cea_mixture_data_ascii( nasa_mixture, Antioch::DefaultFilename::thermo_data() );
  Antioch::CEAEvaluator<Scalar> nasa_thermo( nasa_mixture );

  Antioch::IdealGasMicroThermo< Antioch::NASAEvaluator<Scalar, Antioch::CEACurveFit<Scalar> >,
                                   Scalar > id_thermo( nasa_thermo, chem_mixture );

  // Mass fractions
  std::vector<Scalar> mass_fractions( 5, 0.2 );
  mass_fractions[0] = 0.5;
  mass_fractions[1] = 0.2;
  mass_fractions[2] = 0.1;
  mass_fractions[3] = 0.1;
  mass_fractions[4] = 0.1;

  const Scalar R_N2 = Antioch::Constants::R_universal<Scalar>() / Mm_N2;
  const Scalar R_O2 = Antioch::Constants::R_universal<Scalar>() / Mm_O2;
  const Scalar R_N  = Antioch::Constants::R_universal<Scalar>() / Mm_N;
  const Scalar R_O  = Antioch::Constants::R_universal<Scalar>() / Mm_O;
  const Scalar R_NO = Antioch::Constants::R_universal<Scalar>() / Mm_NO;

  int return_flag = 0;

  Scalar T = 300.1;

  // N2
  return_flag = test_molecule(0,Scalar(2.5),R_N2,nasa_thermo.cv_over_R(Antioch::TempCache<Scalar>(T),0), T, id_thermo, "N2") || return_flag;
  return_flag = test_molecule(1,Scalar(2.5),R_O2,nasa_thermo.cv_over_R(Antioch::TempCache<Scalar>(T),1), T, id_thermo, "O2") || return_flag;
  return_flag = test_molecule(2,Scalar(1.5),R_N, nasa_thermo.cv_over_R(Antioch::TempCache<Scalar>(T),2), T, id_thermo, "N")  || return_flag;
  return_flag = test_molecule(3,Scalar(1.5),R_O, nasa_thermo.cv_over_R(Antioch::TempCache<Scalar>(T),3), T, id_thermo, "O")  || return_flag;
  return_flag = test_molecule(4,Scalar(2.5),R_NO,nasa_thermo.cv_over_R(Antioch::TempCache<Scalar>(T),4), T, id_thermo, "NO") || return_flag;


  return return_flag;
}


int main()
{

  int ierr = (tester<double>()      ||
              tester<long double>() ||
              tester<float>());

  return ierr;
}
