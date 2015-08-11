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

#ifndef ANTIOCH_READ_REACTION_SET_DATA_H
#define ANTIOCH_READ_REACTION_SET_DATA_H

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/units.h"
#include "antioch/parsing_enum.h"

// C++
#include <string>
#include <vector>

namespace Antioch
{

  template <typename NumericType>
  class ASCIIParser;

  template <typename NumericType>
  class ChemKinParser;

  template <typename NumericType>
  class XMLParser;

  template <typename NumericType>
  class ReactionSet;

 /*!\file read_reaction_set_data.h
  *
  * We parse the file here, with an exhaustive
  * unit management. The starting point is the kinetics
  * equation:
  * \f[
  *     \frac{\partial c}{\partial t} = k \prod_{s \in \text{reactants}} c_s^{l_s}
  * \f]
  * with \f$l\f$ the partial order of the reaction with respect to reactants \f$s\f$.
  * We obtain thus
  * \f[
  *     \unit{[k] = [M]^{m - 1} s^{-1}}
  * \f]
  * with \f$m\f$ the order of the reaction. By definition
  * \f[
  *     m = \sum_{s \in \text{reactants}} l_s
  * \f]
  * We are in an elementary processes paradigm, thus for a reactant species \f$s\f$,
  * \f$l_s = -\nu_s\f$ with \f$\nu_s\f$ the stoichiometric coefficient of reactant
  * species \f$s\f$.
  *
  * Example:
  * \f[
  *   \begin{array}{c}
  *      \ce{a A + b B -> c C + d D} \\
  *      m = a + b
  *   \end{array}
  * \f]
  *
  * To this, we consider the kinetics model (they're all included in the
  * Van't Hoff equation):
  * \f[
  *   \alpha(T) = A \left(\frac{T}{\mathrm{T_\text{ref}}}\right)^\beta\exp\left(-\frac{E_a}{\mathrm{R}T} + D T\right)
  * \f]
  *
  * We derive from this all the tests and default units:
  * \f[
  *  \begin{array}{lcccc}\toprule
  *                                & A                                           & \beta & E_a                & D \\\midrule
  *   \text{Elementary}            & \unit{\left(m^3mol^{-1}\right)^{m-1}s^{-1}} & -     & \unit{J\,mol^{-1}} & \unit{K^{-1}} \\
  *   \text{Duplicate}             & \unit{\left(m^3mol^{-1}\right)^{m-1}s^{-1}} & -     & \unit{J\,mol^{-1}} & \unit{K^{-1}} \\
  *   \text{Three body}            & \unit{\left(m^3mol^{-1}\right)^{m}s^{-1}}   & -     & \unit{J\,mol^{-1}} & \unit{K^{-1}} \\
  *   \text{Falloff}\quad k_0      & \unit{\left(m^3mol^{-1}\right)^{m}s^{-1}}   & -     & \unit{J\,mol^{-1}} & \unit{K^{-1}} \\
  *   \text{Falloff}\quad k_\infty & \unit{\left(m^3mol^{-1}\right)^{m-1}s^{-1}} & -     & \unit{J\,mol^{-1}} & \unit{K^{-1}} \\\bottomrule
  *  \end{array}
  * \f]
  * for the Troe falloff, the additionnal parameters are:
  * \f[
  *  \begin{array}{cccc}\toprule
  *    \alpha  & T^*      & T^{**}   & T^{***} \\\midrule
  *     -      & \unit{K} & \unit{K} & \unit{K} \\\bottomrule
  *  \end{array}
  * \f]
  *
  * Thus the reading is made in this fashion:
  *   - read reactants and products, get \f$m\f$
  *   - find default unit of \f$A\f$
  *   - read other parameters
  */
  template<class NumericType>
  void read_reaction_set_data_xml( const std::string& filename,
                                   const bool verbose,
                                   ReactionSet<NumericType>& reaction_set );

  template<class NumericType>
  void read_reaction_set_data_chemkin( const std::string& filename,
                                       const bool verbose,
                                       ReactionSet<NumericType>& reaction_set );


  template<typename NumericType>
  void read_reaction_set_data(const std::string &filename,
                              const bool verbose,
                              ReactionSet<NumericType>& reaction_set,
                              ParsingType type = ASCII );

  template <typename NumericType>
  void verify_unit_of_parameter(Units<NumericType> & default_unit, const std::string & provided_unit,
                                 const std::vector<std::string> & accepted_unit,
                                 const std::string & equation,     const std::string & parameter_name);


  /* ------------------------- Inline Functions -------------------------*/
  template<class NumericType>
  inline
  void read_reaction_set_data_xml( const std::string& filename,
                                   const bool verbose,
                                   ReactionSet<NumericType>& reaction_set )
  {
     read_reaction_set_data<NumericType >(filename,verbose,reaction_set,XML);
  }

  template<class NumericType>
  inline
  void read_reaction_set_data_chemkin( const std::string& filename,
                                       const bool verbose,
                                       ReactionSet<NumericType>& reaction_set )
  {
     read_reaction_set_data<NumericType>(filename,verbose,reaction_set,CHEMKIN);
  }

} // end namespace Antioch

#endif // ANTIOCH_READ_REACTION_SET_DATA_H
