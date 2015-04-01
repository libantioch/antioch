//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014 Paul T. Bauman, Benjamin S. Kirk, Sylvain Plessis,
//                    Roy H. Stonger
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
//--------------------------------------------------------------------------

#ifndef ANTIOCH_PHYSICS_UTILS_H
#define ANTIOCH_PHYSICS_UTILS_H

// Antioch
#include "antioch/physical_set.h"
#include "antioch/transport_mixture.h"

//C++

namespace Antioch
{

  /*! a free fonction to extrapolate to high temperatures
      
     Only the Stockmayer potential imposes a maximal temperature,
     which means only the kinetics theory on transport is
     limited by it.
  */

  // We are explicit on the mixture
  template <typename Physics, typename ThermoEvaluator, typename CoeffType, typename StateType>
  void extrapolate_to_high_temperatures(PhysicalSet<Physics,TransportMixture<ThermoEvaluator,CoeffType> > & physics, const StateType & T_max);


// ---------------------------


  template <typename Physics, typename ThermoEvaluator, typename CoeffType, typename StateType>
  inline
  void extrapolate_to_high_temperatures(PhysicalSet<Physics,TransportMixture<ThermoEvaluator,CoeffType> > & physics, const StateType & T_max)
  {
     extrapolate_T(physics.set(),T_max,typename physical_tag<Physics>::temperature_limitation_type() );
  }

}//end namespace Antioch


#endif //ANTIOCH_PHYSICS_UTILS_H
