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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef ANTIOCH_PHYSICAL_CONSTANTS_H
#define ANTIOCH_PHYSICAL_CONSTANTS_H

namespace Antioch
{
  namespace Constants
  {
    /*!
     * Universal Gas Constant, R, expressed in J/(kmol-K)
     */
    template<typename CoeffType>
    inline
    CoeffType R_universal()
    {
      return 8314.472L;
    }

    /*!
     * Avogadro's number, particles per mole.
     */
    template<typename CoeffType>
    inline
    CoeffType Avogadro()
    {
      return 6.02214179e23L;
    }

  } // end namespace Constants
} // end namespace Antioch

#endif //ANTIOCH_PHYSICAL_CONSTANTS_H
