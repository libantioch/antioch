#ifndef _ELEMENTARY_CHEMICAL_PROCESS_
#define _ELEMENTARY_CHEMICAL_PROCESS_
#include "antioch/Process.hpp"

namespace Antioch{
/*!\file ElementaryProcess.hpp
 * \brief Contains the corresponding class
 *
 * \class ElementaryProcess
 * \brief The simplest chemical reaction
 *
 * It is basically a Process object with
 * a default of Kooij kinetics.
 */

class ElementaryProcess:public Process
{
  public:
/*!\brief Default constructor*/
    ElementaryProcess(const std::string &mod = "Kooij"):
        Process(mod){setKineticsProcess("Elementary process");}
/*!\brief Copy constructor*/
    ElementaryProcess(const ElementaryProcess &rhs):
        Process(*(static_cast<const Process*>(&rhs))){setKineticsProcess("Elementary process");}
/*!\brief Constructing constructor*/
    ElementaryProcess(const std::string &kinMod, const std::vector<ParameterPhy> &pars):
        Process(kinMod,pars){setKineticsProcess("Elementary process");}
/*!\brief Destructor*/
    ~ElementaryProcess(){}

};
}
#endif
