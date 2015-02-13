//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_PARSER_BASE_H
#define ANTIOCH_PARSER_BASE_H

//Antioch

//C++
#include <vector>
#include <string>

namespace Antioch
{
  /*!\class ParserBase

      A parser is an instance related to a file. The parser
      corresponds to a file type (e.g. XML or ChemKin). The
      file HAS to be given in the constructor as the parser
      is indissociable from the file. Set to `true' by default,
      a verbose switch is also available.

      We define here the rule of parsing for
      all parsers. Differences/specificities are described in the
      corresponding files.

      The different things to parse are:
        - the species, name and core characteristics:
                - list of names (std::string) `const std::vector<std::string> species_list() const;'
                - mandatory data (molar mass, heat of formation at 0 K,
                                  number of translational/rotational DOFs,
                                  charge number) `void read_chemical_species(ChemicalMixture<NumericType> & chem_mixture);'
                - vibrational data `void read_vibrational_data(ChemicalMixture<NumericType> & chem_mixture);'
                - electronic data `void read_electronic_data(ChemicalMixture<NumericType> & chem_mixture);'
        - the kinetics: 
                - it requires a boolean to ensure there is a reaction to parse `bool reaction();'
   */

  class ParserBase
  {
        public:
          ParserBase(){}
          virtual ~ParserBase(){}
  };


} // end namespace Antioch


#endif //ANTIOCH_PARSER_BASE_H
