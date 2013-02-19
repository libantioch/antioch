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

// Antioch
#include "antioch/reaction.h"

namespace Antioch
{
  // Forward declarations
  template<typename CoeffType>
  class ReactionSet;

  template<typename CoeffType>
  class ChemicalMixture;
  
  //! Class to handle computing mass source terms for a given ReactionSet.
  /*! This class preallocates work arrays and so *must* be created within a spawned
   *  thread, if running in a threaded environment. It takes a reference to an
   *  already created ReactionSet, so there's little construction penalty.
   */
  template<typename CoeffType=double>
  class Kinetics
  {
  public:

    Kinetics( const ReactionSet<CoeffType>& reaction_set );
    ~Kinetics();

    const ReactionSet<CoeffType>& reaction_set() const;

    //! Compute species production/destruction rates per unit volume in \f$ \left(kg/sec/m^3\right)\f$
    template <typename StateType>
    void compute_mass_sources ( const StateType T,
				const StateType rho,
				const StateType R_mix,
				const std::vector<StateType>& mass_fractions,
				const std::vector<StateType>& molar_densities,
				const std::vector<StateType>& h_RT_minus_s_R,
				std::vector<StateType>& mass_sources );

    unsigned int n_species() const;

    unsigned int n_reactions() const;

  protected:

    const ReactionSet<CoeffType>& _reaction_set;

    const ChemicalMixture<CoeffType>& _chem_mixture;
  };

  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType>
  inline
  const ReactionSet<CoeffType>& Kinetics<CoeffType>::reaction_set() const
  {
    return _reaction_set;
  }

  template<typename CoeffType>
  inline
  unsigned int Kinetics<CoeffType>::n_species() const
  {
    return _chem_mixture.n_species();
  }

  template<typename CoeffType>
  inline
  unsigned int Kinetics<CoeffType>::n_reactions() const
  {
    return _reaction_set.n_reactions();
  }


  template<typename CoeffType>
  inline
  Kinetics<CoeffType>::Kinetics( const ReactionSet<CoeffType>& reaction_set )
    : _reaction_set( reaction_set ),
      _chem_mixture( reaction_set.chemical_mixture() )
  {
    return;
  }


  template<typename CoeffType>
  inline
  Kinetics<CoeffType>::~Kinetics()
  {
    return;
  }


  template<typename CoeffType>
  template<typename StateType>
  inline
  void Kinetics<CoeffType>::compute_mass_sources( const StateType T,
						  const StateType rho,
						  const StateType R_mix,
						  const std::vector<StateType>& mass_fractions,
						  const std::vector<StateType>& molar_densities,
						  const std::vector<StateType>& h_RT_minus_s_R,
						  std::vector<StateType>& mass_sources )
  {
    //! \todo Make these assertions vector-compatible
    // antioch_assert_greater(T, 0.0);
    // antioch_assert_greater(rho, 0.0);
    // antioch_assert_greater(R_mix, 0.0);
    antioch_assert_equal_to( mass_fractions.size(), this->n_species() );
    antioch_assert_equal_to( molar_densities.size(), this->n_species() );
    antioch_assert_equal_to( h_RT_minus_s_R.size(), this->n_species() );
    antioch_assert_equal_to( mass_sources.size(), this->n_species() );

    //! Work arrays for compute mass sources
    // Use a copy constructor here to handle vector StateType sizing
    std::vector<StateType> net_reaction_rates(this->n_reactions(), T);

    std::fill( mass_sources.begin(), mass_sources.end(), 0.0 );
    
    // compute the requisite reaction rates
    this->_reaction_set.compute_reaction_rates( T, rho, R_mix, mass_fractions, molar_densities,
						h_RT_minus_s_R, net_reaction_rates );

    // compute the actual mass sources in kmol/sec/m^3
    for (unsigned int rxn = 0; rxn < this->n_reactions(); rxn++)
      {
	const Reaction<CoeffType>& reaction = this->_reaction_set.reaction(rxn);
	const StateType rate = net_reaction_rates[rxn];
	
	// reactant contributions
	for (unsigned int r = 0; r < reaction.n_reactants(); r++)
	  {
	    const unsigned int r_id = reaction.reactant_id(r);
	    const unsigned int r_stoich = reaction.reactant_stoichiometric_coefficient(r);
	    
	    mass_sources[r_id] -= (static_cast<CoeffType>(r_stoich)*rate);
	  }
	
	// product contributions
	for (unsigned int p=0; p < reaction.n_products(); p++)
	  {
	    const unsigned int p_id = reaction.product_id(p);
	    const unsigned int p_stoich = reaction.product_stoichiometric_coefficient(p);
	    
	    mass_sources[p_id] += (static_cast<CoeffType>(p_stoich)*rate);
	  }
      }

    // finally scale by molar mass
    for (unsigned int s=0; s < this->n_species(); s++)
      {
	mass_sources[s] *= _chem_mixture.M(s);
      }
    
    return;
  }


} // end namespace Antioch

#endif // ANTIOCH_KINETICS_H
