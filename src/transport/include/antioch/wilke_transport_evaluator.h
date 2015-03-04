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

#ifndef ANTIOCH_WILKE_TRANSPORT_EVALUATOR_H
#define ANTIOCH_WILKE_TRANSPORT_EVALUATOR_H

// Antioch
#include "antioch/metaprogramming.h"
#include "antioch/physics_metaprogramming_decl.h"
#include "antioch/temp_cache.h"
#include "antioch/wilke_mixture.h"
#include "antioch/cmath_shims.h"

namespace Antioch
{

  
  template<class Diffusion, class Viscosity, class ThermalConductivity, 
           class Mixture,  // this is a WilkeMixture, but let's stop the templates here...
           class CoeffType=double>
  class WilkeTransportEvaluator
  {
  public:

    WilkeTransportEvaluator( const Mixture& mixture,
                    const Diffusion& diffusion,
                    const Viscosity& viscosity,
                    const ThermalConductivity& conductivity );

    ~WilkeTransportEvaluator();

    //! mixture level diffusion, array of values
    template <typename StateType, typename VectorStateType>
    void D( const StateType& T, const StateType & rho, const StateType & cTot,
            const VectorStateType& mass_fractions, const VectorStateType& mu,
            VectorStateType & ds) const;

    //! mixture level viscosity, one value
    template <typename StateType, typename VectorStateType>
    StateType mu( const StateType& T,
                  const VectorStateType& mass_fractions ) const;

    //! mixture level thermal conduction, one value
    template <typename StateType, typename VectorStateType>
    StateType k( const StateType& T,
                 const VectorStateType& mass_fractions ) const;

    //! mixture level thermal conduction and viscosity, one value
    template <typename StateType, typename VectorStateType>
    void mu_and_k( const StateType& T,
                   const VectorStateType& mass_fractions,
                   StateType& mu, StateType& k ) const;

    //! Thermal conduction and diffusion, helper function
    template <typename StateType, typename VectorStateType>
    void D_and_k(const VectorStateType & mu, const StateType & T, const StateType & rho, 
                 const VectorStateType & mass_fractions, VectorStateType & k, VectorStateType & D) const;

    //! mixture level thermal conduction, viscosity (one value) and diffusion (array of values)
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
    compute_phi( typename Antioch::rebind<VectorStateType,VectorStateType>::type & mu_mu_sqrt,
                 const VectorStateType& chi,
                 const unsigned int s ) const;

    //! Helper function to reduce code duplication.
    /*! Computes the intermediate \mu,\mu matrix needed for Wilke's mixing rule.  */
    template <typename VectorStateType>
    void compute_mu_mu_sqrt( const VectorStateType & mu, 
                             typename Antioch::rebind<VectorStateType,VectorStateType>::type & mu_mu_sqrt) const;

  protected:

    const Mixture& _mixture;

    const Diffusion& _diffusion;

    const Viscosity& _viscosity;

    const ThermalConductivity& _conductivity;

  private:

    WilkeTransportEvaluator();

  };

  template<class Diffusion, class Viscosity, class ThermalConductivity, class Mixture, class CoeffType>
  WilkeTransportEvaluator<Diffusion,Viscosity,ThermalConductivity,Mixture,CoeffType>::WilkeTransportEvaluator( const Mixture& mixture,
                                                                           const Diffusion& diffusion,
                                                                           const Viscosity& viscosity,
                                                                           const ThermalConductivity& conductivity )
    : _mixture(mixture),
      _diffusion(diffusion),
      _viscosity(viscosity),
      _conductivity(conductivity)
  {
    return;
  }

  template<class Diffusion, class Viscosity, class ThermalConductivity, class Mixture, class CoeffType>
  WilkeTransportEvaluator<Diffusion,Viscosity,ThermalConductivity,Mixture,CoeffType>::~WilkeTransportEvaluator()
  {
    return;
  }

  template<class Diffusion, class Viscosity, class ThermalConductivity, class Mixture, class CoeffType>
  template <typename StateType, typename VectorStateType>
  void WilkeTransportEvaluator<Diffusion,Viscosity,ThermalConductivity,Mixture,CoeffType>::D( const StateType& T,
                                                                                        const StateType & rho, const StateType & cTot,
                                                                                        const VectorStateType& mass_fractions,
                                                                                        const VectorStateType& mu,
                                                                                              VectorStateType & ds) const
  {
// first way
// \todo: if not needed, find a way to supress Ds building
// \todo: initialization as full squared matrix, even though half is needed and used
     typename rebind<VectorStateType,VectorStateType>::type Ds(ds.size());
     init_constant(Ds,ds);
     _diffusion.compute_binary_diffusion_matrix(T,cTot,Ds );
     wilke_diffusion_rule(_mixture().transport_mixture().chemical_mixture(), mass_fractions, Ds, ds, typename physical_tag<typename Diffusion::model>::diffusion_species_type ());

// second way
// \todo: if not needed, find a way to supress k building
// \todo: Ds[s][s] not needed here
     StateType k = zero_clone(T);
     for(unsigned int s = 0; s < _mixture.transport_mixture().n_species(); s++)
     {
        _conductivity.compute_thermal_conductivity(s,mu[s],Ds[s][s],T,cTot * _mixture().transport_mixture().chemical_mixture().M(s),k);
        _diffusion.compute_diffusivity(_mixture.thermo_evaluator().cp(TempCache<StateType>(T),s), k, ds[s]); // \todo, solve this by KineticsConditions
     }
  }

  template<class Diffusion, class Viscosity, class ThermalConductivity, class Mixture, class CoeffType>
  template <typename StateType, typename VectorStateType>
  StateType WilkeTransportEvaluator<Diffusion,Viscosity,ThermalConductivity,Mixture,CoeffType>::mu( const StateType& T,
                                                                         const VectorStateType& mass_fractions ) const
  {
    StateType mu_mix = zero_clone(T);
    
    VectorStateType mu  = zero_clone(mass_fractions);
    VectorStateType chi = zero_clone(mass_fractions);

    typename Antioch::rebind<VectorStateType,VectorStateType>::type mu_mu_sqrt(mu.size());
    Antioch::init_constant(mu_mu_sqrt,mu);

    this->compute_mu_chi( T, mass_fractions, mu, chi );

    this->compute_mu_mu_sqrt( mu, mu_mu_sqrt);

    for( unsigned int s = 0; s < _mixture.chem_mixture().n_species(); s++ )
      {
        StateType phi_s = this->compute_phi( mu_mu_sqrt, chi, s );
        
        mu_mix += mu[s]*chi[s]/phi_s;
      }

    return mu_mix;
  }

  template<class Diffusion, class Viscosity, class ThermalConductivity, class Mixture, class CoeffType>
  template <typename StateType, typename VectorStateType>
  StateType WilkeTransportEvaluator<Diffusion,Viscosity,ThermalConductivity,Mixture,CoeffType>::k( const StateType& T,
                                                                        const VectorStateType& mass_fractions ) const
  {
    antioch_assert_equal_to(mass_fractions.size(), _mixture.chem_mixture().n_species());

    StateType k_mix = zero_clone(T);

    VectorStateType mu  = zero_clone(mass_fractions);
    VectorStateType chi = zero_clone(mass_fractions);

    typename Antioch::rebind<VectorStateType,VectorStateType>::type mu_mu_sqrt(mu.size());
    Antioch::init_constant(mu_mu_sqrt,mu);

    this->compute_mu_chi( T, mass_fractions, mu, chi );

    this->compute_mu_mu_sqrt( mu, mu_mu_sqrt);

    for( unsigned int s = 0; s < _mixture.chem_mixture().n_species(); s++ )
      {
        StateType phi_s = this->compute_phi( mu_mu_sqrt, chi, s );
        
        StateType k_s = zero_clone(T);
        _conductivity.compute_thermal_conductivity( s, mu[s], zero_clone(T),T, zero_clone(T), k_s); //\todo, better management

        k_mix += k_s*chi[s]/phi_s;
      }

    return k_mix;
  }

  template<class Diffusion, class Viscosity, class ThermalConductivity, class Mixture, class CoeffType>
  template <typename StateType, typename VectorStateType>
  void WilkeTransportEvaluator<Diffusion,Viscosity,ThermalConductivity,Mixture,CoeffType>::mu_and_k( const StateType& T,
                                                                          const VectorStateType& mass_fractions,
                                                                          StateType& mu_mix,
                                                                          StateType& k_mix ) const
  {
    mu_mix = zero_clone(T);
    k_mix  = zero_clone(T);

    VectorStateType mu  = zero_clone(mass_fractions);
    VectorStateType chi = zero_clone(mass_fractions);

    typename Antioch::rebind<VectorStateType,VectorStateType>::type mu_mu_sqrt(mu.size());
    Antioch::init_constant(mu_mu_sqrt,mu);

    this->compute_mu_chi( T, mass_fractions, mu, chi );

    this->compute_mu_mu_sqrt( mu, mu_mu_sqrt);

    for( unsigned int s = 0; s < _mixture.transport_mixture().n_species(); s++ )
      {
        StateType phi_s = this->compute_phi( mu_mu_sqrt, chi, s );
        
        StateType k_s = zero_clone(T);
        _conductivity.compute_thermal_conductivity( s, mu[s], zero_clone(T),T, zero_clone(T), k_s); //\todo, better management

        mu_mix += mu[s]*chi[s]/phi_s;
        k_mix += k_s*chi[s]/phi_s;
      }

    return;
  }

  template<class Diffusion, class Viscosity, class ThermalConductivity, class Mixture, class CoeffType>
  template <typename StateType, typename VectorStateType>
  void WilkeTransportEvaluator<Diffusion,Viscosity,ThermalConductivity,Mixture,CoeffType>::D_and_k(const VectorStateType & mu, const StateType & T, const StateType & rho, 
                                                                                  const VectorStateType & mass_fractions, VectorStateType & k, VectorStateType & ds) const
  {
     antioch_assert_equal_to(ds.size(),mass_fractions.size());
     antioch_assert_equal_to(ds.size(),_mixture.transport_mixture().n_species());
     antioch_assert_equal_to(mu.size(),_mixture.transport_mixture().n_species());
     const StateType n_molar_mixture = rho / _mixture.transport_mixture().chemical_mixture().M(mass_fractions); // total molar density

// diffusion comes first
// \todo: if not needed, find a way to supress Ds building
// \todo: initialization as full squared matrix, even though half is needed and used
     typename Antioch::rebind<VectorStateType,VectorStateType>::type Ds(mass_fractions.size());
     init_constant(Ds,ds);
     _diffusion.compute_binary_diffusion_matrix(T, n_molar_mixture, Ds);
     wilke_diffusion_rule(_mixture.transport_mixture().chemical_mixture(), mass_fractions, Ds, ds, typename physical_tag<typename Diffusion::model>::diffusion_species_type ());

// thermal conduction
    for(unsigned int s = 0; s < _mixture.transport_mixture().n_species(); s++)
    {
      _diffusion.compute_self_diffusion(s,T,n_molar_mixture,Ds[s][s]);
      _conductivity.compute_thermal_conductivity(s,mu[s],Ds[s][s],T,n_molar_mixture * _mixture.transport_mixture().chemical_mixture().M(s),k[s]);
    }

// diffusion comes last
    for(unsigned int s = 0; s < _mixture.transport_mixture().n_species(); s++)
        _diffusion.compute_diffusivity(rho, _mixture.thermo_evaluator().cp(TempCache<StateType>(T),s), k[s], ds[s] );
  }

  template<class Diffusion, class Viscosity, class ThermalConductivity, class Mixture, class CoeffType>
  template <typename StateType, typename VectorStateType>
  void WilkeTransportEvaluator<Diffusion,Viscosity,ThermalConductivity,Mixture,CoeffType>::mu_and_k_and_D( const StateType& T, const StateType & rho,
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

    typename Antioch::rebind<VectorStateType,VectorStateType>::type mu_mu_sqrt(mu.size());
    Antioch::init_constant(mu_mu_sqrt,mu);

    this->compute_mu_chi( T, mass_fractions, mu, chi );

    this->D_and_k(mu, T, rho , mass_fractions, k, ds );

    this->compute_mu_mu_sqrt( mu, mu_mu_sqrt);

    for( unsigned int s = 0; s < _mixture.transport_mixture().n_species(); s++ )
      {
        StateType phi_s = this->compute_phi( mu_mu_sqrt, chi, s );

        mu_mix += mu[s]*chi[s]/phi_s;
        k_mix += k[s]*chi[s]/phi_s;
      }

    return;
  }

  template<class Diffusion, class Viscosity, class ThermalConductivity, class Mixture, class CoeffType>
  template <typename StateType, typename VectorStateType>
  void WilkeTransportEvaluator<Diffusion,Viscosity,ThermalConductivity,Mixture,CoeffType>::compute_mu_chi( const StateType& T,
                                                                                const VectorStateType& mass_fractions,
                                                                                VectorStateType& mu,
                                                                                VectorStateType& chi ) const
  {
    const StateType M = _viscosity.mixture().M(mass_fractions);

    // Precompute needed quantities
    // chi_s = w_s*M/M_s
    for( unsigned int s = 0; s < _mixture.chem_mixture().n_species(); s++ )
      {
        _viscosity.compute_viscosity(s,T,mu[s]);
        chi[s] = mass_fractions[s]*M/_viscosity.mixture().M(s);
      }

    return;
  }

  template<class Diffusion, class Viscosity, class ThermalConductivity, class Mixture, class CoeffType>
  template <typename VectorStateType>
  typename
  Antioch::value_type<VectorStateType>::type
  WilkeTransportEvaluator<Diffusion,Viscosity,ThermalConductivity,Mixture,CoeffType>::compute_phi(typename Antioch::rebind<VectorStateType,VectorStateType>::type & mu_mu_sqrt,
                                                                        const VectorStateType& chi,
                                                                        const unsigned int s ) const
  {
    using std::sqrt;

    typedef typename
      Antioch::value_type<VectorStateType>::type StateType;

    /* We initialize to the first iterate and loop starting from 1
       since some StateTypes have a hard time initializing from
       a constant. */
    // phi_s = sum_r (chi_r*(1+sqrt(mu_s/mu_r)*(Mr/Ms)^(1/4))^2)/sqrt(8*(1+Ms/Mr))
    const StateType dummy = 1 + mu_mu_sqrt[s][0]*_mixture.Mr_Ms_to_the_one_fourth(0,s);
    StateType phi_s = chi[0]*dummy*dummy/_mixture.denominator(0,s);

    for(unsigned int r = 1; r < _mixture.chem_mixture().n_species(); r++ )
      {
        const StateType numerator = 1 + mu_mu_sqrt[s][r]*_mixture.Mr_Ms_to_the_one_fourth(r,s);
        phi_s += chi[r]*numerator*numerator/_mixture.denominator(r,s);
      }

    return phi_s;
  }

  template<class Diffusion, class Viscosity, class ThermalConductivity, class Mixture, class CoeffType>
  template <typename VectorStateType>
  inline
  void WilkeTransportEvaluator<Diffusion,Viscosity,ThermalConductivity,Mixture,CoeffType>::compute_mu_mu_sqrt( const VectorStateType & mu, 
                                                    typename Antioch::rebind<VectorStateType,VectorStateType>::type & mu_mu_sqrt) const
                                        
  {
    antioch_assert_equal_to(mu.size(),mu_mu_sqrt.size()); // indeed square matrix
    antioch_assert_greater(mu.size(),0); // not empty
#ifndef NDEBUG
    for(unsigned int i = 0; i < mu.size(); i++)
    {
        antioch_assert_equal_to(mu.size(),mu_mu_sqrt[i].size());
    }
#endif

/// first the diagonal
    for(unsigned int s = 0; s < mu.size(); s++)
    {
       mu_mu_sqrt[s][s] = Antioch::constant_clone(mu[s],1);
    }
///then the square roots and the inversion (first row, first column)
   for(unsigned int s = 1; s < mu.size(); s++)
   {
      mu_mu_sqrt[0][s] = ant_sqrt(mu[0]/mu[s]);
      mu_mu_sqrt[s][0] = Antioch::constant_clone(mu[0],1)/mu_mu_sqrt[0][s];
   }
/// now the remaining n-2 matrix, using the fact that mu[i]/mu[j] = mu[i]/mu[0] * mu[0]/mu[j]
  for(unsigned int s = 1; s < mu.size(); s++)
  {
    for(unsigned int l = 1; l < mu.size(); l++)
    {
        if(l == s)continue;
        mu_mu_sqrt[s][l] = mu_mu_sqrt[s][0] * mu_mu_sqrt[0][l];
    }
  }
 }

} // end namespace Antioch

#endif // ANTIOCH_WILKE_TRANSPORT_EVALUATOR_H
