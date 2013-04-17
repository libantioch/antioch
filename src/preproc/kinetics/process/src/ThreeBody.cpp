#include "antioch/ThreeBody.hpp"
#include <stdexcept>

namespace Antioch{
double ThreeBody::getCoefficient(const std::string &molecule) const
{
  try
  {
    return alpha.at(molecule);
  }
  catch(const std::out_of_range& oor)
  {
    return 1.;
  }
  return 1.;
}
}
