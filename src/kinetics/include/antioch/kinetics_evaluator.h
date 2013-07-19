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

#ifndef ANTIOCH_KINETICS_EVALUATOR_H
#define ANTIOCH_KINETICS_EVALUATOR_H

// Antioch
#include "antioch/metaprogramming.h"
#include "antioch/reaction.h"

// C++
#include <vector>

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
  template<typename CoeffType=double, typename StateType=CoeffType>
  class KineticsEvaluator
  {
  public:

    //! Constructor.  Requires a reaction set to be evaluated later,
    //as well as an \p example instantiation of the data type to be
    //used as inputs.  For scalar-valued inputs, "0" is a valid
    //example input.  For vector-valued inputs, an appropriately-sized
    //vector should be used.
    KineticsEvaluator( const ReactionSet<CoeffType>& reaction_set,
                       const StateType& example );

    ~KineticsEvaluator();

    const ReactionSet<CoeffType>& reaction_set() const;

    //! Compute species production/destruction rates per unit volume
    /*! \f$ \left(kg/sec/m^3\right)\f$ */
    template <typename VectorStateType>
    void compute_mass_sources( const StateType& T,
                               const VectorStateType& molar_densities,
                               const VectorStateType& h_RT_minus_s_R,
                               VectorStateType& mass_sources );
    
    //! Compute species production/destruction rate derivatives
    /*! In mass units, e.g. \f$ \frac{\partial \dot{\omega}}{dT}
      [\left(kg/sec/m^3/K\right)]\f$ */
    template <typename VectorStateType>
    void compute_mass_sources_and_derivs( const StateType& T,
                                          const VectorStateType& molar_densities,
                                          const VectorStateType& h_RT_minus_s_R,
                                          const VectorStateType& dh_RT_minus_s_R_dT,
                                          VectorStateType& mass_sources,
                                          VectorStateType& dmass_dT,
                                          std::vector<VectorStateType>& dmass_drho_s );

    //! Compute species molar production/destruction rates per unit volume
    /*! \f$ \left(mole/sec/m^3\right)\f$ */
    template <typename VectorStateType>
    void compute_mole_sources( const StateType& T,
                               const VectorStateType& molar_densities,
                               const VectorStateType& h_RT_minus_s_R,
                               VectorStateType& mole_sources );

    //! Compute species production/destruction rate derivatives
    /*! In mass units, e.g. \f$ \frac{\partial \dot{\omega}}{dT}
      [\left(mole/sec/m^3/K\right)]\f$ */
    template <typename VectorStateType>
    void compute_mole_sources_and_derivs( const StateType& T,
                                          const VectorStateType& molar_densities,
                                          const VectorStateType& h_RT_minus_s_R,
                                          const VectorStateType& dh_RT_minus_s_R_dT,
                                          VectorStateType& mole_sources,
                                          VectorStateType& dmole_dT,
                                          std::vector<VectorStateType>& dmole_dX_s );

    unsigned int n_species() const;

    unsigned int n_reactions() const;

  protected:

    const ReactionSet<CoeffType>& _reaction_set;

    const ChemicalMixture<CoeffType>& _chem_mixture;

    std::vector<StateType> _net_reaction_rates;

    std::vector<StateType> _dnet_rate_dT;

    std::vector<std::vector<StateType> > _dnet_rate_dX_s;
  };

  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType, typename StateType>
  inline
  const ReactionSet<CoeffType>& KineticsEvaluator<CoeffType,StateType>::reaction_set() const
  {
    return _reaction_set;
  }

  template<typename CoeffType, typename StateType>
  inline
  unsigned int KineticsEvaluator<CoeffType,StateType>::n_species() const
  {
    return _chem_mixture.n_species();
  }

  template<typename CoeffType, typename StateType>
  inline
  unsigned int KineticsEvaluator<CoeffType,StateType>::n_reactions() const
  {
    return _reaction_set.n_reactions();
  }


  template<typename CoeffType, typename StateType>
  inline
  KineticsEvaluator<CoeffType,StateType>::KineticsEvaluator
  ( const ReactionSet<CoeffType>& reaction_set,
    const StateType& example )
    : _reaction_set( reaction_set ),
      _chem_mixture( reaction_set.chemical_mixture() ),
      _net_reaction_rates( reaction_set.n_reactions(), example ),
      _dnet_rate_dT( reaction_set.n_reactions(), example ),
      _dnet_rate_dX_s( reaction_set.n_reactions() )
  {

    for( unsigned int r = 0; r < reaction_set.n_reactions(); r++ )
      {
        _dnet_rate_dX_s[r].resize( reaction_set.n_species(), example );
      }

    return;
  }


  template<typename CoeffType, typename StateType>
  inline
  KineticsEvaluator<CoeffType,StateType>::~KineticsEvaluator()
  {
    return;
  }

  template<typename CoeffType, typename StateType>
  template<typename VectorStateType>
  inline
  void KineticsEvaluator<CoeffType,StateType>::compute_mole_sources( const StateType& T,
                                                                     const VectorStateType& molar_densities,
                                                                     const VectorStateType& h_RT_minus_s_R,
                                                                     VectorStateType& mole_sources )
  {
    //! \todo Make these assertions vector-compatible
    // antioch_assert_greater(T, 0.0);
    antioch_assert_equal_to( molar_densities.size(), this->n_species() );
    antioch_assert_equal_to( h_RT_minus_s_R.size(), this->n_species() );
    antioch_assert_equal_to( mole_sources.size(), this->n_species() );

    /*! \todo Do we need to really initialize this? */
    Antioch::set_zero(_net_reaction_rates);

    Antioch::set_zero(mole_sources);

    // compute the requisite reaction rates
    this->_reaction_set.compute_reaction_rates( T, molar_densities,
                                                h_RT_minus_s_R, _net_reaction_rates );

    // compute the actual mole sources in kmol/sec/m^3
    for (unsigned int rxn = 0; rxn < this->n_reactions(); rxn++)
      {
        const Reaction<CoeffType>& reaction = this->_reaction_set.reaction(rxn);
        const StateType rate = _net_reaction_rates[rxn];
        
        // reactant contributions
        for (unsigned int r = 0; r < reaction.n_reactants(); r++)
          {
            const unsigned int r_id = reaction.reactant_id(r);
            const unsigned int r_stoich = reaction.reactant_stoichiometric_coefficient(r);
            
            mole_sources[r_id] -= (static_cast<CoeffType>(r_stoich)*rate);
          }
        
        // product contributions
        for (unsigned int p=0; p < reaction.n_products(); p++)
          {
            const unsigned int p_id = reaction.product_id(p);
            const unsigned int p_stoich = reaction.product_stoichiometric_coefficient(p);
            
            mole_sources[p_id] += (static_cast<CoeffType>(p_stoich)*rate);
          }
      }

    return;
  }


  template<typename CoeffType, typename StateType>
  template<typename VectorStateType>
  inline
  void KineticsEvaluator<CoeffType,StateType>::compute_mass_sources( const StateType& T,
                                                                     const VectorStateType& molar_densities,
                                                                     const VectorStateType& h_RT_minus_s_R,
                                                                     VectorStateType& mass_sources )
  {
    // Quantities asserted in compute_mole_sources call
    this->compute_mole_sources( T, molar_densities, h_RT_minus_s_R, mass_sources );

    // finally scale by molar mass
    for (unsigned int s=0; s < this->n_species(); s++)
      {
        mass_sources[s] *= _chem_mixture.M(s);
      }
    
    return;
  }

  template<typename CoeffType, typename StateType>
  template<typename VectorStateType>
  inline
  void KineticsEvaluator<CoeffType,StateType>::compute_mole_sources_and_derivs( const StateType& T,
                                                                                const VectorStateType& molar_densities,
                                                                                const VectorStateType& h_RT_minus_s_R,
                                                                                const VectorStateType& dh_RT_minus_s_R_dT,
                                                                                VectorStateType& mole_sources,
                                                                                VectorStateType& dmole_dT,
                                                                                std::vector<VectorStateType>& dmole_dX_s )
  {
    //! \todo Make these assertions vector-compatible
    // antioch_assert_greater(T, 0.0);
    antioch_assert_equal_to( molar_densities.size(), this->n_species() );
    antioch_assert_equal_to( h_RT_minus_s_R.size(), this->n_species() );
    antioch_assert_equal_to( dh_RT_minus_s_R_dT.size(), this->n_species() );
    antioch_assert_equal_to( mole_sources.size(), this->n_species() );
    antioch_assert_equal_to( dmole_dT.size(), this->n_species() );
    antioch_assert_equal_to( dmole_dX_s.size(), this->n_species() );
#ifdef NDEBUG
#else
    for (unsigned int s=0; s < this->n_species(); s++)
      {
        antioch_assert_equal_to( dmole_dX_s[s].size(), this->n_species() );
      }
#endif
    
    /*! \todo Do we need to really initialize these? */
    Antioch::set_zero(_net_reaction_rates);
    Antioch::set_zero(_dnet_rate_dT);

    Antioch::set_zero(mole_sources);
    Antioch::set_zero(dmole_dT);
    for (unsigned int s=0; s < this->n_species(); s++)
      {
        Antioch::set_zero(dmole_dX_s[s]);

        /*! \todo Do we need to really initialize this? */
        Antioch::set_zero(_dnet_rate_dX_s[s]);
      }

    // compute the requisite reaction rates
    this->_reaction_set.compute_reaction_rates_and_derivs( T, molar_densities, 
                                                           h_RT_minus_s_R, dh_RT_minus_s_R_dT,
                                                           _net_reaction_rates,
                                                           _dnet_rate_dT, 
                                                           _dnet_rate_dX_s );
    // compute the actual mole sources in kmol/sec/m^3
    for (unsigned int rxn = 0; rxn < this->n_reactions(); rxn++)
      {

        const Reaction<CoeffType>& reaction = this->_reaction_set.reaction(rxn);

        /*! \todo Are these going to get optimized out? Should we remove them? */
        const StateType rate = _net_reaction_rates[rxn];
        const StateType drate_dT = _dnet_rate_dT[rxn];
        const VectorStateType drate_dX_s = _dnet_rate_dX_s[rxn];

        // reactant contributions
        for (unsigned int r = 0; r < reaction.n_reactants(); r++)
          {
            const unsigned int r_id = reaction.reactant_id(r);
            const unsigned int r_stoich = reaction.reactant_stoichiometric_coefficient(r);
            
            mole_sources[r_id] -= (static_cast<CoeffType>(r_stoich)*rate);

            // d/dT rate contributions
            dmole_dT[r_id] -= (static_cast<CoeffType>(r_stoich)*drate_dT);

            // d(.m)/dX_s rate contributions
            for (unsigned int s=0; s < this->n_species(); s++)
              {
                dmole_dX_s[r_id][s] -= (static_cast<CoeffType>(r_stoich)*drate_dX_s[s]);
              }
          }
        
        // product contributions
        for (unsigned int p=0; p < reaction.n_products(); p++)
          {
            const unsigned int p_id = reaction.product_id(p);
            const unsigned int p_stoich = reaction.product_stoichiometric_coefficient(p);
            
            mole_sources[p_id] += (static_cast<CoeffType>(p_stoich)*rate);

            // d/dT rate contributions
            dmole_dT[p_id] += (static_cast<CoeffType>(p_stoich)*drate_dT);

            // d/dX_s rate contributions
            for (unsigned int s=0; s < this->n_species(); s++)
              {
                dmole_dX_s[p_id][s] += (static_cast<CoeffType>(p_stoich)*drate_dX_s[s]);
              }
          }
      }

    return;
  }


  template<typename CoeffType, typename StateType>
  template <typename VectorStateType>
  inline
  void KineticsEvaluator<CoeffType,StateType>::compute_mass_sources_and_derivs( const StateType& T,
                                                                                const VectorStateType& molar_densities,
                                                                                const VectorStateType& h_RT_minus_s_R,
                                                                                const VectorStateType& dh_RT_minus_s_R_dT,
                                                                                VectorStateType& mass_sources,
                                                                                VectorStateType& dmass_dT,
                                                                                std::vector<VectorStateType>& dmass_drho_s )
  {
    // Asserts are in compute_mole_sources
    this->compute_mole_sources_and_derivs( T, molar_densities, h_RT_minus_s_R, dh_RT_minus_s_R_dT,
                                           mass_sources, dmass_dT, dmass_drho_s );
    
    // Convert from mole units to mass units
    for (unsigned int s=0; s < this->n_species(); s++)
      {
        mass_sources[s] *= _chem_mixture.M(s);
        dmass_dT[s] *= _chem_mixture.M(s);

        for (unsigned int t=0; t < this->n_species(); t++)
          {
            dmass_drho_s[s][t] *= _chem_mixture.M(s)/_chem_mixture.M(t);
          }

      }

    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_KINETICS_EVALUATOR_H
