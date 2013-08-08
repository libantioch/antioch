//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
#ifndef _PREPROC_PHYSICAL_CONSTANTS_
#define _PREPROC_PHYSICAL_CONSTANTS_

#include "antioch/ParameterPhy.hpp"

namespace Antioch{
/* \namespace PhyCon
 * \brief Useful physical constants.
 *
 * This namespace defines useful physical constants.
 */

namespace PhyCon{
/*!\brief double for kinetics reference temperature, in K*/
static const double Tkin(300.);
/*!\brief double for ideal gas constant, in J/mol/K*/
static const double Rgas(8.3144621);
/*!\brief Ideal gas contant in J/mol/K from NIST 2010*/
static const ParameterPhy R("R",PhyCon::Rgas,0.0000075,CORE_UNCERTAINTY_TYPE_ABSOLUTE,"J/mol/K");
/*!\brief Tref for Kooij equation, set to 300 K*/
static const ParameterPhy TrefKin("Tref",PhyCon::Tkin,0.,CORE_UNCERTAINTY_TYPE_NONE,"K");
/*!\brief Useful for unit management*/
static const ParameterPhy NAvogadro("Avogadro constant",6.02214129e23,2.7e16,CORE_UNCERTAINTY_TYPE_ABSOLUTE,"mol-1");
};
}
#endif
