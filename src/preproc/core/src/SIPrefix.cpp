//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
#include "antioch/SIPrefix.hpp"


/************************
 * SIPrefixes methods   *
 ************************/

namespace Antioch{
SIPrefixes &SIPrefixes::operator=(const SIPrefixes &rhs)
{
  if(this == &rhs){return *this;}
  _symbol = rhs.symbol();
  _value=rhs.value();
  return *this;
}
}
