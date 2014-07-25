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
 *
 * This library is intended for performances, thus the reference temperature is taken equal
 * to one in order to skip the division step (\f$\frac{T}{\mathrm{T_0}}\f$). This
 * is \e not an option, this is an \e obligation: the kinetics equations are coded without the
 * division step.
 */
namespace Antioch
{
  namespace KineticsModel
  {

    enum KineticsModel { CONSTANT = 0,   // A
                         HERCOURT_ESSEN, // A * T^beta
                         BERTHELOT,      // A * exp(D*T)
                         ARRHENIUS,      // A * exp(-Ea/T)
                         BHE,            // A * T^beta * exp(D*T)
                         KOOIJ,          // A * T^beta * exp(-Ea/T)
                         VANTHOFF,       // A * T^beta * exp(-Ea/T + D*T)
                         PHOTOCHEM };    // int_0^\infty f(\lambda)\sigma(\lambda) d\lambda = const(T)

    template<typename CoeffType>
    CoeffType Tref()
    {
      return 1.0; // this HAS to stay this way because it is hard-coded for performances (see eq. above)
    }

  } // end namespace KineticsModel

} // end namespace Antioch

#endif // ANTIOCH_REACTION_ENUM_H
