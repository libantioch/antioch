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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef ANTIOCH_WILKE_EVALUATOR_H
#define ANTIOCH_WILKE_EVALUATOR_H

// Antioch
#include "antioch/metaprogramming.h"
#include "antioch/wilke_mixture.h"
#include "antioch/mixture_viscosity.h"
#include "antioch/wilke_transport_evaluator.h"
#include "antioch/physics_placeholder.h"
#include "antioch/antioch_asserts.h"

namespace Antioch
{

  
  template<class Viscosity, class ThermalConductivity, class CoeffType=double>
  class WilkeEvaluator:
        public WilkeTransportEvaluator<PhysicsPlaceholder, Viscosity, ThermalConductivity, WilkeMixture<CoeffType>,CoeffType>
  {
  public:

    WilkeEvaluator( const WilkeMixture<CoeffType>& mixture,
                    const Viscosity& viscosity,
                    const ThermalConductivity& conductivity );

    ~WilkeEvaluator();

  private:

    WilkeEvaluator();

  };

  template<class Viscosity, class ThermalConductivity, class CoeffType>
  WilkeEvaluator<Viscosity,ThermalConductivity,CoeffType>::WilkeEvaluator( const WilkeMixture<CoeffType>& mixture,
                                                                           const Viscosity& viscosity,
                                                                           const ThermalConductivity& conductivity )
    : WilkeTransportEvaluator<PhysicsPlaceholder,Viscosity,ThermalConductivity,WilkeMixture<CoeffType>,CoeffType>(mixture,PhysicsPlaceholder(),viscosity,conductivity)
  {
    antioch_deprecated();
    return;
  }

  template<class Viscosity, class ThermalConductivity, class CoeffType>
  WilkeEvaluator<Viscosity,ThermalConductivity,CoeffType>::~WilkeEvaluator()
  {
    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_WILKE_EVALUATOR_H
