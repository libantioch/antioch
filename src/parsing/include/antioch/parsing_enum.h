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

namespace Antioch
{
  enum ParsingKey{REACTION_DATA = 0,
                  REACTION,
                  REVERSIBLE,
                  ID,
                  EQUATION,
                  CHEMICAL_PROCESS,
                  KINETICS_MODEL,
                  REACTANTS,
                  PRODUCTS,
                  PREEXP,
                  POWER,
                  ACTIVATION_ENERGY,
                  BERTHELOT_COEFFICIENT,
                  TREF,
                  HV_LAMBDA,
                  HV_CROSS_SECTION,
                  UNIT,
                  EFFICIENCY,
                  TROE_FALLOFF,
                  TROE_F_ALPHA,
                  TROE_F_TS,
                  TROE_F_TSS,
                  TROE_F_TSSS
                 };
}
