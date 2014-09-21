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

#ifndef ANTIOCH_TRANSPORT_SPECIES_H
#define ANTIOCH_TRANSPORT_SPECIES_H

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/physical_constants.h"
#include "antioch/Lennard_Jones_potential.h"

// C++
#include <vector>
#include <iostream>
#include <limits>

namespace Antioch
{

  typedef unsigned int Species;

  //! Class to encapsulate data relevant for transport for each chemical species
  /*!
   * This class is designed to store information relevant to the transport of a chemical species.
   * All the data stored is constant for each species, so we store const for each
   * variable. The idea is that this will be placed inside TransportMixture, which
   * will be a singleton, dependant of ChemicalMixture.
   */
  template<typename CoeffType=double>
  class TransportSpecies
  {
  public:
    
    //! Constructor
    TransportSpecies(const Species   name, 
                     const CoeffType LJ_depth,
                     const CoeffType LJ_diameter,
                     const CoeffType dipole_moment,
                     const CoeffType polarizability,
                     const CoeffType rot_relax,
                     const CoeffType mass);

    //! Default constructor.
    /*!
     * This is technically required for any
     * std::map value type (or operator[] breaks, at least).  But,
     * we never actually want to create a SpeciesChemistry
     * implicitly, so we throw an error if this is ever used.
     */      
    TransportSpecies();

    //! Destructor
    ~TransportSpecies();

    //! returns a descriptive name for this species.
    const Species species() const;

    //!returns the Lennard-Jones depth in (K), this is \f$\frac{\epsilon}{\mathrm{k_B}}\f$.
    CoeffType LJ_depth() const;

    //!returns the Lennard-Jones diameter in (Angström).
    CoeffType LJ_diameter() const;

    //! return Lennard-Jones potential
    LennardJonesPotential<CoeffType> & LJ() const;

    //! returns dipole moment in units of [D]
    CoeffType dipole_moment() const;

    //! Returns polarizability in units of [Angström^3]
    CoeffType polarizability() const;

    //! Returns rotational relaxation collision number at 298 K, no unit
    CoeffType rotational_relaxation() const;

    //! Returns molecular mass in kg
    CoeffType M() const;

    //! boolean testing polarity
    bool polar() const;
    
    //! Formatted print
    /*!
     * Defaults to \p std::cout.
     */
    void print(std::ostream &os = std::cout) const;

    //!Formatted print 
    /*!
     * Allows you to do std::cout << object << std::endl;
     */
    friend std::ostream& operator<<( std::ostream& os,
                                     const TransportSpecies<CoeffType>& species )
    {
      species.print(os);
      return os;
    }

  protected:

    //! chemical species
    const Species _name;

    //!Lennard-Jones potential
    LennardJonesPotential<CoeffType> _LJ;

    //! Dipole moment in units of [D]
    CoeffType _dipole_moment;

    //! Polarizability in units of [Angström^3]
    CoeffType _polarizability;

    //! Rotational relaxation
    CoeffType _rotational_relaxation;

    //! molar mass in kg/mol
    CoeffType _mass;

  }; // class TransportSpecies


  /* ------------------------- Friend Functions ------------------------- */

  

  /* ------------------------- Inline Functions ------------------------- */

  template<typename CoeffType>
  inline
  const Species TransportSpecies<CoeffType>::species() const 
  { return _name; }

  template<typename CoeffType>
  inline
  CoeffType TransportSpecies<CoeffType>::LJ_depth() const
  { return _LJ.depth(); }

  template<typename CoeffType>
  inline
  CoeffType TransportSpecies<CoeffType>::LJ_diameter() const 
  { return _LJ.diameter(); }

  template<typename CoeffType>
  inline
  LennardJonesPotential<CoeffType> & TransportSpecies<CoeffType>::LJ() const
  {
    return _LJ;
  }

  template<typename CoeffType>
  inline
  CoeffType TransportSpecies<CoeffType>::dipole_moment() const 
  { return _dipole_moment; }

  template<typename CoeffType>
  inline
  CoeffType TransportSpecies<CoeffType>::polarizability() const 
  { return _polarizability; }

  template<typename CoeffType>
  inline
  CoeffType TransportSpecies<CoeffType>::rotational_relaxation() const 
  { return _rotational_relaxation; }

  template<typename CoeffType>
  inline
  CoeffType TransportSpecies<CoeffType>::M() const 
  { return _mass; }

  template<typename CoeffType>
  inline
  bool TransportSpecies<CoeffType>::polar() const 
  {
     return (_dipole_moment > std::numeric_limits<CoeffType>::epsilon());
  }

  template<typename CoeffType>
  inline
  TransportSpecies<CoeffType>::TransportSpecies( const Species name,
                                                 const CoeffType LJ_depth,
                                                 const CoeffType LJ_diameter,
                                                 const CoeffType dipole_moment,
                                                 const CoeffType polarizability,
                                                 const CoeffType rotational_relaxation,
                                                 const CoeffType mass)
    : _name                  (name),
      _LJ                    (LJ_depth,LJ_diameter),
      _dipole_moment         (dipole_moment),
      _polarizability        (polarizability),
      _rotational_relaxation (rotational_relaxation),
      _mass                  (mass)
  {
    return;
  }


  template<typename CoeffType>
  inline
  TransportSpecies<CoeffType>::TransportSpecies()
    : _name(-1), 
      _LJ(0.,0.),
      _polarizability(0.),
      _rotational_relaxation(0.)
  {
    antioch_error();
    return;
  }


  template<typename CoeffType>
  inline
  TransportSpecies<CoeffType>::~TransportSpecies()
  {
    return;
  }


  template<typename CoeffType>
  inline
  void TransportSpecies<CoeffType>::print (std::ostream &os) const
  {
    os << " -----------------------------\n"
       << "| Species enum " << this->species() << '\n'
       << "| Transport data "             << '\n'
       << " -----------------------------\n"
       << std::scientific
       << "  LJ depth       = " << this->LJ_depth()              << '\n'
       << "  LJ diameter    = " << this->LJ_diameter()           << '\n'
       << "  polarizability = " << this->polarizability()        << '\n'
       << "  rot relax      = " << this->rotational_relaxation() << '\n';

    os << '\n';
    
    return;
  }


} // end namespace Antioch

#endif // ANTIOCH_CHEMICAL_SPECIES_H
