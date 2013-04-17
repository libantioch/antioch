//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#ifndef _ANTIOCH_KINETICS_FACTORY_H
#define _ANTIOCH_KINETICS_FACTORY_H

//Antioch
#include "antioch/kinetics_enum.h"
#include "antioch/kinetics_type.h"
#include "antioch/hercourtessen_rate.h"
#include "antioch/berthelot_rate.h"
#include "antioch/arrhenius_rate.h"
#include "antioch/berthelothercourtessen_rate.h"
#include "antioch/kooij_rate.h"
#include "antioch/vanthoff_rate.h"

namespace Antioch{

KineticsType *kinetics(const KinMod::KinMod kin)
{
  switch(kin)
  {
     case(KinMod::HERCOURT_ESSEN):
     {
     }
        break;
     case(KinMod::BERTHELOT):
     {
     }
        break;
     case(KinMod::ARRHENIUS):
     {
     }
        break;
     case(KinMod::BHE):
     {
     }
        break;
     case(KinMod::KOOIJ):
     {
     }
        break;
     case(KinMod::VANTHOFF):
     {
     }
        break;
  }
}

}

#endif
