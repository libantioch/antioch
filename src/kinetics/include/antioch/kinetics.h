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

#ifndef ANTIOCH_KINETICS_H
#define ANTIOCH_KINETICS_H

// C++
#include <vector>

namespace Antioch
{
  // Forward declarations
  template<class NumericType>
  class ReactionSet;

  template<class NumericType>
  class ChemicalMixture;
  
  //! Class to handle computing mass source terms for a given ReactionSet.
  /*! This class preallocates work arrays and so *must* be created within a spawned
   *  thread, if running in a threaded environment. It takes a reference to an
   *  already created ReactionSet, so there's little construction penalty.
   */
  template<class NumericType>
  class Kinetics
  {
  public:

    Kinetics( const ReactionSet<NumericType>& reaction_set );
    ~Kinetics();

    const ReactionSet<NumericType>& reaction_set() const;

    //! Compute species production/destruction rates per unit volume in \f$ \left(kg/sec/m^3\right)\f$
    void compute_mass_sources ( const NumericType T,
				const NumericType rho,
				const NumericType R_mix,
				const std::vector<NumericType>& mass_fractions,
				const std::vector<NumericType>& molar_densities,
				const std::vector<NumericType>& h_RT_minus_s_R,
				std::vector<NumericType>& mass_sources );

    unsigned int n_species() const;

    unsigned int n_reactions() const;

  protected:

    const ReactionSet<NumericType>& _reaction_set;

    const ChemicalMixture<NumericType>& _chem_mixture;

    //! Work arrays for compute mass sources
    std::vector<NumericType> _net_reaction_rates;
  };

  /* ------------------------- Inline Functions -------------------------*/
  template<class NumericType>
  inline
  const ReactionSet<NumericType>& Kinetics<NumericType>::reaction_set() const
  {
    return _reaction_set;
  }

  template<class NumericType>
  inline
  unsigned int Kinetics<NumericType>::n_species() const
  {
    return _chem_mixture.n_species();
  }

  template<class NumericType>
  inline
  unsigned int Kinetics<NumericType>::n_reactions() const
  {
    return _reaction_set.n_reactions();
  }

} // end namespace Antioch

#endif // ANTIOCH_KINETICS_H
