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
#ifndef _GENERAL_CHEMICAL_METHOD_
#define _GENERAL_CHEMICAL_METHOD_

#include <vector>
#include <string>
/*!\file GeneralChemistry.hpp
 * \brief General public functions useful for chemistry
 *
 * This file contains public methods to perform general
 * operations, as splitting a molecule into atoms and
 * isomere comparisons. Two levels of isomere is provided (warning, these
 * levels are not chemically relevant, they are purely computing relevant, see
 * concerned functions for details):
 *   - strict isomere, see bool isStrictIsomere(const std::string&, const std::string&),
 *   - ``loose'' isomere, see bool isIsomere(const std::string&, const std::string&).
 */

/*!\brief Split the given molecule into atoms
 *
 * There is no check that the splitting give actual
 * real atoms. An atom is characterized by its number A, the 
 * energy state and its number of occurences in
 * the molecule. The considered format is [A]Xx(Es) with
 * A the number of nucleons, Xx the atom's name and Es
 * its energy state, such as 4S or 2D. If the
 * first characters are figures, they are considered to
 * be the number of nucleons of the first atom.
 */
std::vector<std::string> splitInAtoms(const std::string &molecule);

/*!\brief Molecule comparator, strict version
 *
 * Two isomeres are two molecules with the same atoms. Yet
 * there are several levels of isomere that can be considered,
 * as the structure of two isomeres can differ widely.
 *
 * Thus a strict isomere will be two isomeres (two molecules
 * with the same atoms) which formulae provide the atoms
 * in the same order, \e e.g. \e CH2 and CHH. CH2 and HCH
 * will not be strict isomeres.
 */
bool isStrictIsomere(const std::string &mol1, const std::string &mol2);

/*!\brief Molecule comparator
 *
 * Two isomeres are molecules composed of the same atoms, regardless
 * of the structure.
 */
bool isIsomere(const std::string &mol1, const std::string &mol2);

#endif
