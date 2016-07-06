//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014-2016 Paul T. Bauman, Benjamin S. Kirk,
//                         Sylvain Plessis, Roy H. Stonger
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

#ifndef ANTIOCH_CONDUCTIVITY_TRAITS_H
#define ANTIOCH_CONDUCTIVITY_TRAITS_H

#include "antioch/eucken_thermal_conductivity.h"
#include "antioch/kinetics_theory_thermal_conductivity.h"

namespace Antioch
{
  template<typename ConductivityModel>
  struct ConductivityTraits;

  template<typename ThermoEvaluator>
  struct ConductivityTraits<EuckenThermalConductivity<ThermoEvaluator> >
  {
    static bool const requires_diffusion = false;
  };

  template <typename ThermoEvaluator, typename CoeffType>
  struct ConductivityTraits<KineticsTheoryThermalConductivity<ThermoEvaluator,CoeffType> >
  {
    static bool const requires_diffusion = true;
  };

} // end namespace Antioch

#endif // ANTIOCH_CONDUCTIVITY_TRAITS_H
