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

#ifndef ANTIOCH_PHYSICAL_CONSTANTS_H
#define ANTIOCH_PHYSICAL_CONSTANTS_H

//Antioch
#include "antioch/units.h"

namespace Antioch
{
  namespace Constants
  {
    /*!
     * Universal Gas Constant, R, expressed in J/(mol-K)
     */
    template<typename CoeffType>
    inline
    CoeffType R_universal()
    {
      return 8.3144621L;
    }

    /*!
     * Avogadro's number, particles per mole.
     */
    template<typename CoeffType>
    inline
    CoeffType Avogadro()
    {
      return 6.02214129e23L;
    }

    /*!
     * Universal Gas Constant unit
     */
    template<typename CoeffType>
    inline
    Units<CoeffType> R_universal_unit()
    {
        return Units<CoeffType>("J/mol/K");
    }

    /*!
     * Avogadro's number unit
     */
    template<typename CoeffType>
    inline
    Units<CoeffType> Avogadro_unit()
    {
        return Units<CoeffType>("mol-1");
    }

    /*
     * Planck's constant, m2.kg/s (J.s)
     */
    template<typename CoeffType>
    inline
    CoeffType Planck_constant()
    {
      return 6.62606957e-34;
    }

    /*!
     * Planck's constant unit
     */
    template<typename CoeffType>
    inline
    Units<CoeffType> Planck_constant_unit()
    {
        return Units<CoeffType>("m2.kg/s");
    }

    /*!
     * light celerity, m/s
     */
    template<typename CoeffType>
    inline
    CoeffType light_celerity()
    {
      return 2.99792458e8;
    }

    /*!
     * light celerity unit
     */
    template<typename CoeffType>
    inline
    Units<CoeffType> light_celerity_unit()
    {
        return Units<CoeffType>("m/s");
    }

    /*!
     * Boltzmann constant
     * 1.380 6488 x 10-23 J/K
     * (http://physics.nist.gov/cgi-bin/cuu/Value?k)
     */
    template<typename CoeffType>
    inline
    CoeffType Boltzmann_constant()
    {
       return 1.3806488e-23;
    }

    /*!
     * Boltzmann constant unit
     */
    template<typename CoeffType>
    inline
    Units<CoeffType> Boltzmann_constant_unit()
    {
        return Units<CoeffType>("J/K");
    }

    /*!
     * vacuum permeability
     *  4 * pi * 10-7 m.kg.s2.A2
     */
    template<typename CoeffType>
    inline
    CoeffType vacuum_permeability()
    {
       return 4.L * pi<CoeffType>() * 1e-7;
    }

    /*!
     * vacuum permeability
     */
    template<typename CoeffType>
    inline
    Units<CoeffType> vacuum_permeability_unit()
    {
        return Units<CoeffType>("m.kg.s2.A2");
    }

    /*!
     * vacuum permittivity
     * 1/(vacuum_permeability * light_celerity^2) A2.s4/kg/m3
     */
    template<typename CoeffType>
    inline
    CoeffType vacuum_permittivity()
    {
       return 1.L/ (vacuum_permeability<CoeffType>() * light_celerity<CoeffType>() 
                                                     * light_celerity<CoeffType>());
    }

    /*!
     * vacuum permeability unit
     */
    template<typename CoeffType>
    inline
    Units<CoeffType> vacuum_permittivity_unit()
    {
        return Units<CoeffType>("A2.s4/kg/m3");
    }

  } // end namespace Constants
} // end namespace Antioch

#endif //ANTIOCH_PHYSICAL_CONSTANTS_H
