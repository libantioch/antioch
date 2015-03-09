//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#include "antioch/chemical_mixture.h"
#include "antioch/transport_mixture.h"
#include "antioch/ascii_parser.h"
#include "antioch/ideal_gas_micro_thermo.h"
#include "antioch/cea_curve_fit.h"
#include "antioch/nasa_mixture.h"
#include "antioch/nasa_evaluator.h"
#include "antioch/thermo_handler.h"

#include <iomanip>
#include <string>

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
    Antioch::NASAThermoMixture<Scalar, Antioch::NASA9CurveFit<Scalar> > nasa_mixture( chem_mix);
    Antioch::NASAEvaluator<Scalar, Antioch::NASA9CurveFit<Scalar> >     nasa_thermo( nasa_mixture );

    typedef Antioch::NASAEvaluator<Scalar, Antioch::NASA9CurveFit<Scalar> > NASAThermoType;

//// micro
    Antioch::IdealGasMicroThermo<NASAThermoType,Scalar>  micro_thermo(nasa_thermo,chem_mix);

    typedef Antioch::IdealGasMicroThermo<NASAThermoType,Scalar>  StatThermoType;

//// eval
    Antioch::ThermoHandler<Scalar,NASAThermoType,StatThermoType> thermo(nasa_thermo,micro_thermo);

    typedef Antioch::ThermoHandler<Scalar,NASAThermoType,StatThermoType> ThermoEval;

//transport
///// mixture and stat thermo for conduction
        // read
    Antioch::ASCIIParser<Scalar> parser(filename,true);

    parser.set_ignored_columns(std::vector<unsigned int>(1,0));

    Antioch::TransportMixture<ThermoEval,Scalar> tran_mixture( chem_mix, thermo, static_cast<Antioch::ParserBase<Scalar>*>(&parser));
  
    return tran_mixture.transport_species().empty(); // false is winner
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
