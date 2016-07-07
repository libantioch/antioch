//antioch
#include <antioch/vector_utils.h>
#include <antioch/antioch_asserts.h>
#include <antioch/chemical_species.h>
#include <antioch/chemical_mixture.h>
#include <antioch/reaction_set.h>
#include <antioch/chemkin_parser.h>
#include <antioch/xml_parser.h>
#include "antioch/nasa_mixture.h"
#include "antioch/nasa_mixture_parsing.h"
#include <antioch/read_reaction_set_data.h>
#include <antioch/kinetics_evaluator.h>
#include <antioch/nasa_evaluator.h>

int main()
{

   int n_species = 8;
   std::string reaction_filename = "five_rxn.xml";
   std::string thermo_filename = "nasa7_thermo.xml";

   /**********************************
    ***
    ***   set up Antioch to do chemistry
    ***
    ***********************************/

    std::vector<double> Y(n_species, 0.0);
    double T = 450.0;
 
    // Define species 
    std::vector<std::string> species_str_list;
    species_str_list.reserve(n_species);
    species_str_list.push_back("H2");
    species_str_list.push_back("O2");
    species_str_list.push_back("H");
    species_str_list.push_back("O");
    species_str_list.push_back("OH");
    species_str_list.push_back("HO2");
    species_str_list.push_back("H2O");
    species_str_list.push_back("N2");
       
    // Get chemistry for species involved in this reaction
    Antioch::ChemicalMixture<double> chem_mixture( species_str_list );

    // Get thermodynamic data
    Antioch::NASAThermoMixture<double, Antioch::NASA7CurveFit<double> > nasa_mixture( chem_mixture );
    Antioch::read_nasa_mixture_data( nasa_mixture, thermo_filename, Antioch::XML );

    // Prepare for chemical reactions
    Antioch::ReactionSet    <double> reaction_set( chem_mixture );
    Antioch::read_reaction_set_data_xml<double>( reaction_filename, true, reaction_set );

    // Set up reactions and thermodynamics 
    Antioch::KineticsEvaluator<double> kinetics(reaction_set, 0); // Reactions
    Antioch::NASAEvaluator<double, Antioch::NASA7CurveFit<double> > thermo(nasa_mixture); // Thermodynamics

    // Prepare Antioch to get chemistry and thermodynamic information
    typedef Antioch::TempCache<double> Cache;

    // Get some thermodynamic data for each species
    std::vector<double> h_RT_minus_s_R(n_species);
    thermo.h_RT_minus_s_R(Cache(T),h_RT_minus_s_R);

    // Perform kinetics calculations 
    // (basically just gives RHS of ODEs for each species)
    // (Note that temperature is passed in but the RHS for the
    //  temperature equation is NOT computed by Antioch.  Just
    //  need temperature in reaction rate calculations.)
    std::vector<double> mole_sources(n_species,0.0);
    std::vector<double> molar_densities(n_species,0.0);
    kinetics.compute_mole_sources(
        T,               // this is temperature (needed for calculations)
        molar_densities, // current concentrations (in moles)
        h_RT_minus_s_R,  // exponent in equilibrium constant for reverse rates
        mole_sources);   // RHS

    // Energy equation computations (thermo)
    std::vector<double> h(n_species, 0.);  // Enthalpy for each species
    std::vector<double> cv(n_species, 0.); // Specific volume for each species
    std::vector<double> cp(n_species, 0.); // Specific pressure for each species

    double Qnum = 0.0;
    double Qden = 0.0;

    for (int s = 0; s < n_species; s++)
    { // Get numerator and denominator in energy equation (sum of species)
        // Get enthalpy and convert to molar from mass basis
        // Note that R_universal/Rs = Ws where Ws is the molecular weight of species s
        h[s] = thermo.h(Cache(T),s)/chem_mixture.R(s)* Antioch::Constants::R_universal<double>();
        // get cp and convert from mass to molar
        cp[s] = thermo.cp(Cache(T),s)/chem_mixture.R(s)* Antioch::Constants::R_universal<double>();
        // Numerator in energy equation: u_s * dx_s/dt
        Qnum += h[s] * mole_sources[s];
        // Denominator in energy equation
        Qden += cp[s] * molar_densities[s];
    }

  return 0;
}
