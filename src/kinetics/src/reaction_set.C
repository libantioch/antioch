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

// This class
#include "antioch/reaction_set.h"

// Antioch
#include "antioch/chemical_species.h"

namespace Antioch
{
  template<class NumericType>
  ReactionSet<NumericType>::ReactionSet( const ChemicalMixture<NumericType>& chem_mixture )
    : _chem_mixture(chem_mixture)
  {
    return;
  }

  template<class NumericType>
  ReactionSet<NumericType>::~ReactionSet()
  {
    return;
  }
  
  template<class NumericType>
  void ReactionSet<NumericType>::compute_reaction_rates ( const NumericType T,
							  const NumericType rho,
							  const NumericType R_mix,
							  const std::vector<NumericType>& mass_fractions,
							  const std::vector<NumericType>& molar_densities,
							  const std::vector<NumericType>& h_RT_minus_s_R,
							  std::vector<NumericType>& net_reaction_rates ) const
  {
    antioch_assert_equal_to( net_reaction_rates.size(), this->n_reactions() );
    antioch_assert_greater(T, 0.0);
    antioch_assert_greater(rho, 0.0);
    antioch_assert_greater(R_mix, 0.0);
    antioch_assert_equal_to( mass_fractions.size(), this->n_species() );
    antioch_assert_equal_to( molar_densities.size(), this->n_species() );
    antioch_assert_equal_to( h_RT_minus_s_R.size(), this->n_species() );

    // useful constants
    const NumericType P0    = 1.e5; // standard pressure
    const NumericType RT    = R_mix*T;
    const NumericType P0_RT = P0 / RT; // used to transform equilibrium constant from pressure units

    // compute reaction forward rates & other reaction-sized arrays
    for (unsigned int rxn=0; rxn<this->n_reactions(); rxn++)
      {
	const Reaction<NumericType>& reaction = this->reaction(rxn);

	NumericType kfwd = (reaction.forward_rate())(T);

	NumericType keq = reaction.equilibrium_constant( P0_RT, h_RT_minus_s_R );

	const NumericType kbkwd = kfwd/keq;

	net_reaction_rates[rxn] = reaction.compute_rate_of_progress( molar_densities, kfwd, kbkwd );
      }
    
    return;
  }

  template<class NumericType>
  void ReactionSet<NumericType>::print(std::ostream& os) const
  {
    os << "# Number of reactions: " << this->n_reactions() << "\n";

    for (unsigned int r=0; r < this->n_reactions(); r++)
      {
	os << "# " << r << '\n'
	   << this->reaction(r) << "\n";
      }

    return;
  }
  
  /* ------------------------- Instantiate ------------------------- */
  template class ReactionSet<double>;

} // end namespace Antioch
