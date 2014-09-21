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

#ifndef ANTIOCH_PURE_SPECIES_THERMAL_CONDUCTIVITY_BUILDING_H
#define ANTIOCH_PURE_SPECIES_THERMAL_CONDUCTIVITY_BUILDING_H

// Antioch
#include "antioch/pure_species_thermal_conductivity.h"
#include "antioch/tranport_mixture.h"
#include "antioch/physical_set.h"

// C++
#include <iostream>
#include <vector>


namespace Antioch
{
  template<class ThermoEvaluator, class NumericType>
  void build_pure_species_thermal_conductivity( PhysicalSet<PureSpeciesThermalConductivity<ThermoEvaluator,NumericType> , TransportMixture<NumericType> >& k)

// ----------------------------------------- //

  template<class NumericType>
  void build_pure_species_thermal_conductivity( PhysicalSet<PureSpeciesThermalConductivity<ThermoEvaluator,NumericType>, TransportMixture<NumericType> >& k)
  {
       for(unsigned int s = 0; s < k.mixture().n_species(); s++)
       {
           Initializer<PureSpeciesThermalConductivity<ThermoEvaluator,NumericType> > init;
           init.t        = k.mixture().thermo();
           init.Z_298K   = k.mixture().transport_species()[s].rotational_relaxation();
           init.M        = k.mixture().transport_species()[s].M();
           init.LJ_depth = k.mixture().transport_species()[s].LJ_depth();
           k.add_model(k.mixture().species_inverse_name_map().at(s),init);
       }
  }

}

#endif
