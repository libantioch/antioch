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

#ifndef ANTIOCH_PARSING_ENUM_H
#define ANTIOCH_PARSING_ENUM_H

namespace Antioch
{
  enum ParsingType{ASCII = 0,
                   XML,
                   CHEMKIN};

  enum ParsingKey{SPECIES_SET = 0,
                  SPECIES_DATA,
                  SPECIES,
                  THERMO,
                  PHASE_BLOCK,
//
                  TMIN,
                  TMAX,
                  NASADATA,
                  NASA7,
                  NASA9,
//
                  REACTION_DATA,
                  REACTION,
                  REVERSIBLE,
                  ID,
                  EQUATION,
                  CHEMICAL_PROCESS,
                  KINETICS_MODEL,
                  REACTANTS,
                  PRODUCTS,
                  FORWARD_ORDER,
                  BACKWARD_ORDER,
                  PREEXP,
                  POWER,
                  ACTIVATION_ENERGY,
                  BERTHELOT_COEFFICIENT,
                  TREF,
                  HV_LAMBDA,
                  HV_CROSS_SECTION,
                  UNIT,
                  EFFICIENCY,
                  FALLOFF_LOW,
                  FALLOFF_LOW_NAME,
                  TROE_FALLOFF,
                  TROE_F_ALPHA,
                  TROE_F_TS,
                  TROE_F_TSS,
                  TROE_F_TSSS
                 };

// GRI30 compatibility
  enum GRI30Comp{FALLOFF = 0,
                 TROE
                };

  enum ParsingUnit{
                    MOL_WEIGHT = 0,
                    MASS_ENTHALPY
                   };

}

#endif
