#include "antioch/LindemannFalloff.hpp"

namespace Antioch{
ParameterPhy LindemannFalloff::FLind("Lindemann falloff F factor",1.0,0.0,CORE_UNCERTAINTY_TYPE_NONE,"");

LindemannFalloff::LindemannFalloff(const LindemannFalloff &rhs)
{
  setKineticsProcess("Lindemann falloff");
  init(rhs.getKineticsModel(),rhs.getParametersFromRate(rates[0]),rhs.getParametersFromRate(rates[1]));
  setTemperature(rhs.getTemperature());
  setConcentration(rhs.getConcentration());
}
}
