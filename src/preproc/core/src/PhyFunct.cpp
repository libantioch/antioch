//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
#include "antioch/PhyFunct.hpp"

namespace Antioch
{

void PhyFunct::showAll(std::ostream &out) const
{
  out << "#PhyFunct" << std::endl;
  out << "name: " << name << std::endl;
  out << "Abscissa" << std::endl;
  x.showAll(out);
  out << "Ordinate" << std::endl;
  y.showAll(out);
}

}
