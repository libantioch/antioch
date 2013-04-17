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

#ifndef ANTIOCH_KINETICS_ENUM_H
#define ANTIOCH_KINETICS_ENUM_H

/*!
 * The parameters are reduced parameters. The complete Van't Hoff
 * equation is (all the others are deductible from it)
 * \f$k(T) = A * \left(\frac{T}{\mathrm{T_0}}\right)^\beta \exp\left(-\frac{E_a}{\mathrm{R}T} + D T\right)\f$
 * with \f$\mathrm{T_0}\f$ being a reference temperature and \f$\mathrm{R}\f$
 * the ideal gas constant.
 * The reduced parameters are:
 * \f[
 *     \begin{array}{ccc}\toprule
 *     \text{Parameter} & \text{reduced parameter} \\\midrule
 *     A         &  A \mathrm{T_0}^\beta \\
 *     \beta     &  \beta      \\
 *     E_a       &  \frac{E_a}{\mathrm{R}}\\
 *     D         &  D\\\bottomrule
 *     \end{array}
 * \f]
 */
namespace Antioch
{
  namespace KinMod
  {

    enum KinMod { HERCOURT_ESSEN = 0,// A * T^beta
			BERTHELOT,   // A * exp(D*T)
                        ARRHENIUS,   // A * exp(-Ea/T)
                        BHE,         // A * T^beta * exp(D*T)
                        KOOIJ,       // A * T^beta * exp(-Ea/T)
                        VANTHOFF };  // A * T^beta * exp(-Ea/T + D*T)

  } // end namespace KinMod
} // end namespace Antioch

#endif // ANTIOCH_REACTION_ENUM_H
