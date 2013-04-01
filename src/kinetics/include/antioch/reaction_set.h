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

#ifndef ANTIOCH_REACTION_SET_H
#define ANTIOCH_REACTION_SET_H

// Antioch
#include "antioch/chemical_mixture.h"
#include "antioch/reaction.h"

// C++
#include <iostream>
#include <vector>

namespace Antioch
{
  /*!
   * This class encapsulates all the reaction mechanisms considered in a
   * chemical nonequilibrium simulation.
   */
  template<typename CoeffType=double>
  class ReactionSet
  {

  public:

    //! Constructor.
    ReactionSet( const ChemicalMixture<CoeffType>& chem_mixture );

    ~ReactionSet();

    //! \returns the number of species.
    unsigned int n_species() const;
     
    //! \returns the number of reactions.
    unsigned int n_reactions() const;
     
    //! Add a reaction to the system.
    void add_reaction(const Reaction<CoeffType>& reaction);

    //! \returns a constant reference to reaction \p r.
    const Reaction<CoeffType>& reaction(const unsigned int r) const;

    const ChemicalMixture<CoeffType>& chemical_mixture() const;

    //! Compute the rates of progress for each reaction
    template <typename StateType, typename VectorStateType, typename VectorReactionsType>
    void compute_reaction_rates( const StateType& T,
				 const StateType& rho,
				 const StateType& R_mix,
				 const VectorStateType& mass_fractions,
				 const VectorStateType& molar_densities,
				 const VectorStateType& h_RT_minus_s_R,
				 VectorReactionsType& net_reaction_rates ) const;

    //! Formatted print, by default to \p std::cout.
    void print( std::ostream& os = std::cout ) const;
     
    //! Formatted print.
    friend std::ostream& operator<<( std::ostream& os, const ReactionSet<CoeffType>& rset )
    {
      rset.print(os);
      return os;
    }

  private:
     
    ReactionSet();
     
    const ChemicalMixture<CoeffType>& _chem_mixture;

    std::vector<Reaction<CoeffType> > _reactions;

  };
  
  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType>
  inline
  unsigned int ReactionSet<CoeffType>::n_species() const
  {
    return _chem_mixture.n_species();
  }

  template<typename CoeffType>
  inline
  unsigned int ReactionSet<CoeffType>::n_reactions() const
  {
    return _reactions.size();
  }

  template<typename CoeffType>
  inline
  void ReactionSet<CoeffType>::add_reaction(const Reaction<CoeffType>& reaction)
  {
    _reactions.push_back(reaction);
    
    // and make sure it is initialized!
    _reactions.back().initialize();

    return;
  }
  
  template<typename CoeffType>
  inline
  const Reaction<CoeffType>& ReactionSet<CoeffType>::reaction(const unsigned int r) const      
  {
    antioch_assert_less(r, this->n_reactions());
    return _reactions[r];
  }

  template<typename CoeffType>
  inline
  const ChemicalMixture<CoeffType>& ReactionSet<CoeffType>::chemical_mixture() const
  {
    return _chem_mixture;
  }


  template<typename CoeffType>
  inline
  ReactionSet<CoeffType>::ReactionSet( const ChemicalMixture<CoeffType>& chem_mixture )
    : _chem_mixture(chem_mixture)
  {
    return;
  }


  template<typename CoeffType>
  inline
  ReactionSet<CoeffType>::~ReactionSet()
  {
    return;
  }
  

  template<typename CoeffType>
  template<typename StateType, typename VectorStateType, typename VectorReactionsType>
  inline
  void ReactionSet<CoeffType>::compute_reaction_rates ( const StateType& T,
							const StateType& rho,
							const StateType& R_mix,
							const VectorStateType& mass_fractions,
							const VectorStateType& molar_densities,
							const VectorStateType& h_RT_minus_s_R,
							VectorReactionsType& net_reaction_rates ) const
  {
    antioch_assert_equal_to( net_reaction_rates.size(), this->n_reactions() );

    //!\todo Make these assertions vector-compatible
    // antioch_assert_greater(T, 0.0);
    // antioch_assert_greater(rho, 0.0);
    // antioch_assert_greater(R_mix, 0.0);

    antioch_assert_equal_to( mass_fractions.size(), this->n_species() );
    antioch_assert_equal_to( molar_densities.size(), this->n_species() );
    antioch_assert_equal_to( h_RT_minus_s_R.size(), this->n_species() );

    // useful constants
    const CoeffType P0    = 1.e5; // standard pressure
    const StateType RT    = R_mix*T;
    const StateType P0_RT = P0 / RT; // used to transform equilibrium constant from pressure units

    // compute reaction forward rates & other reaction-sized arrays
    for (unsigned int rxn=0; rxn<this->n_reactions(); rxn++)
      {
	const Reaction<CoeffType>& reaction = this->reaction(rxn);

	StateType kfwd = (reaction.forward_rate())(T);

	StateType keq = reaction.equilibrium_constant( P0_RT, h_RT_minus_s_R );

	const StateType kbkwd = kfwd/keq;

	net_reaction_rates[rxn] = reaction.compute_rate_of_progress( molar_densities, kfwd, kbkwd );
      }
    
    return;
  }


  template<typename CoeffType>
  inline
  void ReactionSet<CoeffType>::print(std::ostream& os) const
  {
    os << "# Number of reactions: " << this->n_reactions() << "\n";

    for (unsigned int r=0; r < this->n_reactions(); r++)
      {
	os << "# " << r << '\n'
	   << this->reaction(r) << "\n";
      }

    return;
  }
 
} // end namespace Antioch

#endif // ANTIOCH_REACTION_SET_H
