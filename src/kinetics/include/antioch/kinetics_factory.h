//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
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
