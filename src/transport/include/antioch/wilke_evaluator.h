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

#ifndef ANTIOCH_WILKE_EVALUATOR_H
#define ANTIOCH_WILKE_EVALUATOR_H

// Antioch
#include "antioch/metaprogramming.h"
#include "antioch/wilke_mixture.h"
#include "antioch/mixture_viscosity.h"

namespace Antioch
{

  
  template<class Mixture>
  class WilkeEvaluator
  {
  public:

    WilkeEvaluator( const Mixture & mixture);

    ~WilkeEvaluator();

    template <typename StateType, typename VectorStateType>
    StateType mu( const StateType& T,
                  const VectorStateType& mass_fractions ) const;

    template <typename StateType, typename VectorStateType>
    StateType k( const StateType& T,
                 const VectorStateType& mass_fractions, const StateType &rho = 0) const;

    template <typename StateType, typename VectorStateType>
    void D( const StateType& T, const StateType & rho,
            const VectorStateType& mass_fractions, VectorStateType & Ds ) const;

    template <typename StateType, typename VectorStateType>
    void mu_and_k_and_D( const StateType& T, const StateType & rho,
                   const VectorStateType& mass_fractions,
                   StateType& mu, StateType& k, VectorStateType & Ds ) const;

    //! Helper function to reduce code duplication.
    /*! Populates species viscosities and the intermediate \chi variable
        needed for Wilke's mixing rule. */
    template <typename StateType, typename VectorStateType>
    void compute_mu_chi( const StateType& T,
                         const VectorStateType& mass_fractions,
                         VectorStateType& mu,
                         VectorStateType& chi ) const;

    //! Helper function to reduce code duplication.
    /*! Computes the intermediate \phi variable needed for Wilke's mixing rule.  */
    template <typename VectorStateType>
    typename
    Antioch::value_type<VectorStateType>::type
    compute_phi( const VectorStateType& mu,
                 const VectorStateType& chi,
                 const unsigned int s ) const;

  protected:

    const Mixture & _mixture;

  private:

    WilkeEvaluator();

  };

  template<typename Mixture>
  WilkeEvaluator<Mixture>::WilkeEvaluator( const Mixture & mixture)
    : _mixture(mixture)
  {
    return;
  }

  template<typename Mixture>
  WilkeEvaluator<Mixture>::~WilkeEvaluator()
  {
    return;
  }

  template<typename Mixture>
  template <typename StateType, typename VectorStateType>
  StateType WilkeEvaluator<Mixture>::mu( const StateType& T,
                                         const VectorStateType& mass_fractions ) const
  {
    StateType mu_mix = zero_clone(T);
    
    VectorStateType mu  = zero_clone(mass_fractions);
    VectorStateType chi = zero_clone(mass_fractions);
    
    this->compute_mu_chi( T, mass_fractions, mu, chi );

    for( unsigned int s = 0; s < _mixture.transport_mixture().n_species(); s++ )
      {
        StateType phi_s = this->compute_phi( mu, chi, s );
        
        mu_mix += mu[s]*chi[s]/phi_s;
      }

    return mu_mix;
  }

  template<typename Mixture>
  template <typename StateType, typename VectorStateType>
  StateType WilkeEvaluator<Mixture>::k( const StateType& T, 
                                        const VectorStateType& mass_fractions, const StateType & rho ) const
  {
    antioch_assert_equal_to(mass_fractions.size(), _mixture.transport_mixture().n_species());

    StateType k_mix = zero_clone(T);

    VectorStateType mu  = zero_clone(mass_fractions);
    VectorStateType chi = zero_clone(mass_fractions);

    this->compute_mu_chi( T, mass_fractions, mu, chi );

    const StateType n_molar_mixture = rho / _mixture.transport_mixture().chemical_mixture().M(mass_fractions); // total molar density

    for( unsigned int s = 0; s < _mixture.transport_mixture().n_species(); s++ )
      {
        StateType phi_s = this->compute_phi( mu, chi, s );
        
        StateType k_s = _mixture.k( s, mu[s], T, n_molar_mixture );

        k_mix += k_s*chi[s]/phi_s;
      }

    return k_mix;
  }

  template<typename Mixture>
  template <typename StateType, typename VectorStateType>
  void WilkeEvaluator<Mixture>::mu_and_k_and_D( const StateType& T, const StateType & rho,
                                                const VectorStateType& mass_fractions,
                                                StateType& mu_mix,
                                                StateType& k_mix,
                                                VectorStateType & ds ) const
  {
    mu_mix = zero_clone(T);
    k_mix  = zero_clone(T);
    ds = zero_clone(mass_fractions);

    VectorStateType mu  = zero_clone(mass_fractions);
    VectorStateType k  = zero_clone(mass_fractions);
    VectorStateType chi = zero_clone(mass_fractions);

    this->compute_mu_chi( T, mass_fractions, mu, chi );

    _mixture.D_and_k(mu, T, rho , mass_fractions, k, ds );

    for( unsigned int s = 0; s < _mixture.transport_mixture().n_species(); s++ )
      {
        StateType phi_s = this->compute_phi( mu, chi, s );

        mu_mix += mu[s]*chi[s]/phi_s;
        k_mix += k[s]*chi[s]/phi_s;
      }

    return;
  }

  template <typename Mixture>
  template <typename StateType, typename VectorStateType>
  inline
  void WilkeEvaluator<Mixture>::D( const StateType& T, const StateType & rho, const VectorStateType& mass_fractions, VectorStateType & Ds ) const
  {
      VectorStateType mu = zero_clone(mass_fractions);
      _mixture.mu(T,mu);

      _mixture.D(T, rho, rho / _mixture.transport_mixture().chemical_mixture().M(mass_fractions), mass_fractions, mu, Ds);
  }

  template<typename Mixture>
  template <typename StateType, typename VectorStateType>
  void WilkeEvaluator<Mixture>::compute_mu_chi( const StateType& T,
                                                const VectorStateType& mass_fractions,
                                                VectorStateType& mu,
                                                VectorStateType& chi ) const
  {
    const StateType M = _mixture.transport_mixture().chemical_mixture().M(mass_fractions);

    // Precompute needed quantities
    // chi_s = w_s*M/M_s
    for( unsigned int s = 0; s < _mixture.transport_mixture().n_species(); s++ )
      {
        mu[s] = _mixture.mu(s,T);
        chi[s] = mass_fractions[s]*M/_mixture.transport_mixture().chemical_mixture().M(s);
      }

    return;
  }

  template<typename Mixture>
  template <typename VectorStateType>
  typename
  Antioch::value_type<VectorStateType>::type
  WilkeEvaluator<Mixture>::compute_phi( const VectorStateType& mu,
                                        const VectorStateType& chi,
                                        const unsigned int s ) const
  {

    typedef typename
      Antioch::value_type<VectorStateType>::type StateType;

    /* We initialize to the first iterate and loop starting from 1
       since some StateTypes have a hard time initializing from
       a constant. */
    // phi_s = sum_r (chi_r*(1+sqrt(mu_s/mu_r)*(Mr/Ms)^(1/4))^2)/sqrt(8*(1+Ms/Mr))
    const StateType dummy = 1 + ant_sqrt(mu[s]/mu[0])*_mixture.Mr_Ms_to_the_one_fourth(0,s);
    StateType phi_s = chi[0]*dummy*dummy/_mixture.denominator(0,s);

    for(unsigned int r = 1; r < _mixture.transport_mixture().n_species(); r++ )
      {
        const StateType numerator = 1 + ant_sqrt(mu[s]/mu[r])*_mixture.Mr_Ms_to_the_one_fourth(r,s);
        phi_s += chi[r]*numerator*numerator/_mixture.denominator(r,s);
      }

    return phi_s;
  }

} // end namespace Antioch

#endif // ANTIOCH_WILKE_EVALUATOR_H
