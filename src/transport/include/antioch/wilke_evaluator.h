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
#include "antioch/metaprogramming_decl.h"
#include "antioch/wilke_mixture.h"
#include "antioch/mixture_viscosity.h"

namespace Antioch
{

  /*!
   * WilkeEvaluator, small object to easily access
   * transport calculations.
   * KC stands for KineticsConditions or StateType,
   * it adapts at the preprocessing time.
   *
   */
  
  template<class Mixture>
  class WilkeEvaluator
  {
  public:

    WilkeEvaluator( const Mixture & mixture);

    ~WilkeEvaluator();

    template <typename KC, typename VectorStateType>
    typename value_type<VectorStateType>::type 
              mu( const KC& cond,
                  const VectorStateType& mass_fractions ) const;

    template <typename KC, typename VectorStateType>
    typename value_type<VectorStateType>::type 
              k( const KC& cond,
                 const VectorStateType& mass_fractions, 
                 const typename value_type<VectorStateType>::type &rho = 0) const;

    template <typename KC, typename VectorStateType>
    void D( const KC& cond, const typename value_type<VectorStateType>::type & rho,
            const VectorStateType& mass_fractions, VectorStateType & Ds ) const;

    template <typename StateType, typename VectorStateType, typename KC>
    void mu_and_k_and_D( const KC& cond, const StateType & rho,
                   const VectorStateType& mass_fractions,
                   StateType& mu, StateType& k, VectorStateType & Ds ) const;

    //! Helper function to reduce code duplication.
    /*! Populates species viscosities and the intermediate \chi variable
        needed for Wilke's mixing rule. */
    template <typename KC, typename VectorStateType>
    void compute_mu_chi( const KC& cond,
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
  inline
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
  template <typename KC, typename VectorStateType>
  inline
  typename value_type<VectorStateType>::type 
            WilkeEvaluator<Mixture>::mu( const KC& cond,
                                         const VectorStateType& mass_fractions ) const
  {

      // convenient
    typedef typename value_type<VectorStateType>::type StateType;

    typename constructor_or_reference<const KineticsConditions<StateType,VectorStateType>, const KC>::type  //either (KineticsConditions<> &) or (KineticsConditions<>)
                conditions(cond);

    StateType mu_mix = zero_clone(conditions.T());
    
    VectorStateType mu  = zero_clone(mass_fractions);
    VectorStateType chi = zero_clone(mass_fractions);
    
    this->compute_mu_chi( conditions, mass_fractions, mu, chi );

    for( unsigned int s = 0; s < _mixture.transport_mixture().n_species(); s++ )
      {
        StateType phi_s = this->compute_phi( mu, chi, s );
        
        mu_mix += mu[s]*chi[s]/phi_s;
      }

    return mu_mix;
  }

  template<typename Mixture>
  template <typename KC, typename VectorStateType>
  typename value_type<VectorStateType>::type 
            WilkeEvaluator<Mixture>::k( const KC& cond, 
                                        const VectorStateType& mass_fractions, 
                                        const typename value_type<VectorStateType>::type & rho ) const
  {
    antioch_assert_equal_to(mass_fractions.size(), _mixture.transport_mixture().n_species());

    typedef typename value_type<VectorStateType>::type StateType;

    typename constructor_or_reference<const KineticsConditions<StateType,VectorStateType>, const KC>::type  //either (KineticsConditions<> &) or (KineticsConditions<>)
                conditions(cond);

    StateType k_mix = zero_clone(conditions.T());

    VectorStateType mu  = zero_clone(mass_fractions);
    VectorStateType chi = zero_clone(mass_fractions);

    this->compute_mu_chi( conditions, mass_fractions, mu, chi );

    const StateType n_molar_mixture = rho / _mixture.transport_mixture().chemical_mixture().M(mass_fractions); // total molar density

    for( unsigned int s = 0; s < _mixture.transport_mixture().n_species(); s++ )
      {
        StateType phi_s = this->compute_phi( mu, chi, s );
        
        StateType k_s = _mixture.k( s, mu[s], conditions, n_molar_mixture );

        k_mix += k_s*chi[s]/phi_s;
      }

    return k_mix;
  }

  template<typename Mixture>
  template <typename StateType, typename VectorStateType, typename KC>
  void WilkeEvaluator<Mixture>::mu_and_k_and_D( const KC& cond, const StateType & rho,
                                                const VectorStateType& mass_fractions,
                                                StateType& mu_mix,
                                                StateType& k_mix,
                                                VectorStateType & ds ) const
  {

    typename constructor_or_reference<const KineticsConditions<StateType,VectorStateType>, const KC>::type  //either (KineticsConditions<> &) or (KineticsConditions<>)
                conditions(cond);

     
    /*! \todo Do we need to really initialize this? */
    set_zero(mu_mix);
    set_zero(k_mix);
    set_zero(ds);

    VectorStateType mu  = zero_clone(mass_fractions);
    VectorStateType k   = zero_clone(mass_fractions);
    VectorStateType chi = zero_clone(mass_fractions);

    this->compute_mu_chi( conditions, mass_fractions, mu, chi );

    _mixture.D_and_k(mu, conditions, rho , mass_fractions, k, ds );

    for( unsigned int s = 0; s < _mixture.transport_mixture().n_species(); s++ )
      {
        StateType phi_s = this->compute_phi( mu, chi, s );

        mu_mix += mu[s]*chi[s]/phi_s;
        k_mix += k[s]*chi[s]/phi_s;
      }

    return;
  }

  template <typename Mixture>
  template <typename KC, typename VectorStateType>
  inline
  void WilkeEvaluator<Mixture>::D( const KC& cond, 
                                   const typename value_type<VectorStateType>::type & rho, 
                                   const VectorStateType& mass_fractions, VectorStateType & Ds ) const
  {
    typedef typename value_type<VectorStateType>::type StateType;

    typename constructor_or_reference<const KineticsConditions<StateType,VectorStateType>, const KC>::type  //either (KineticsConditions<> &) or (KineticsConditions<>)
                conditions(cond);

      VectorStateType mu = zero_clone(mass_fractions);
      _mixture.mu(conditions,mu);

      _mixture.D(conditions, rho, rho / _mixture.transport_mixture().chemical_mixture().M(mass_fractions), mass_fractions, mu, Ds);
  }

  template<typename Mixture>
  template <typename KC, typename VectorStateType>
  void WilkeEvaluator<Mixture>::compute_mu_chi( const KC& cond,
                                                const VectorStateType& mass_fractions,
                                                VectorStateType& mu,
                                                VectorStateType& chi ) const
  {
    typedef typename value_type<VectorStateType>::type StateType;

    typename constructor_or_reference<const KineticsConditions<StateType,VectorStateType>, const KC>::type  //either (KineticsConditions<> &) or (KineticsConditions<>)
                conditions(cond);

    const StateType M = _mixture.transport_mixture().chemical_mixture().M(mass_fractions);

    // Precompute needed quantities
    // chi_s = w_s*M/M_s
    for( unsigned int s = 0; s < _mixture.transport_mixture().n_species(); s++ )
      {
        mu[s] = _mixture.mu(s,conditions);
        chi[s] = mass_fractions[s] * M / _mixture.transport_mixture().chemical_mixture().M(s);
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
