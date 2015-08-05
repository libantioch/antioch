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
#include "antioch/mixture_averaged_transport_evaluator.h"
#include "antioch/antioch_asserts.h"
#include "antioch/eucken_thermal_conductivity_building.h"

namespace Antioch
{
  //! Deprecated. Use WilkeTransportEvaluator instead.
  template<class MixtureViscosity, class ThermalConductivity, class CoeffType=double>
  class WilkeEvaluator
  {
  public:

    WilkeEvaluator( const WilkeMixture<CoeffType>& mixture,
                    const MixtureViscosity& viscosity,
                    const ThermalConductivity& conductivity );

    ~WilkeEvaluator();


    template <typename StateType, typename VectorStateType>
    StateType mu( const StateType& T,
                  const VectorStateType& mass_fractions ) const {return _wilke_eval->mu(T,mass_fractions);}

    template <typename StateType, typename VectorStateType>
    StateType k( const StateType& T,
                 const VectorStateType& mass_fractions ) const {return _wilke_eval->k(T,mass_fractions);}

    template <typename StateType, typename VectorStateType>
    void mu_and_k( const StateType& T,
                   const VectorStateType& mass_fractions,
                   StateType& mu, StateType& k ) const {_wilke_eval->mu_and_k(T,mass_fractions,mu,k);}

    //! Helper function to reduce code duplication.
    /*! Populates species viscosities and the intermediate \chi variable
        needed for Wilke's mixing rule. */
    template <typename StateType, typename VectorStateType>
    void compute_mu_chi( const StateType& T,
                         const VectorStateType& mass_fractions,
                               VectorStateType& mu,
                               VectorStateType& chi ) const {_wilke_eval->compute_mu_chi(T,mass_fractions,mu,chi);}

     //! Helper function to reduce code duplication.
     /*! Computes the intermediate \phi variable needed for Wilke's mixing rule. */
     template <typename VectorStateType>
     typename
       Antioch::value_type<VectorStateType>::type
        compute_phi( const VectorStateType& mu,
                     const VectorStateType& chi,
                     const unsigned int s ) const;


  private:

    WilkeEvaluator();

    TransportMixture<CoeffType>* _transport_mixture;

    WilkeTransportMixture<CoeffType>* _wilke_mixture;

    //! This is dummy.
    /*! WilkeEvaluator doesn't support diffusion, but the new WilkeTransportEvaluator does
        so we create a dummy diffusion class, don't bother building it and pass it along
        to WilkeTransportEvaluator. */
    MixtureDiffusion<ConstantLewisDiffusivity<CoeffType>,CoeffType>* _diffusion;

    MixtureConductivity<ThermalConductivity,CoeffType>* _conductivity;

    typename ThermalConductivity::micro_thermo_type* _micro_thermo;

    typedef  WilkeTransportEvaluator<ConstantLewisDiffusivity<CoeffType>,
                                     typename MixtureViscosity::species_viscosity_type,
                                     ThermalConductivity,
                                     CoeffType> Evaluator;

    Evaluator* _wilke_eval;

  };

  template<class Viscosity, class ThermalConductivity, class CoeffType>
  WilkeEvaluator<Viscosity,ThermalConductivity,CoeffType>::WilkeEvaluator( const WilkeMixture<CoeffType>& mixture,
                                                                           const Viscosity& viscosity,
                                                                           const ThermalConductivity& /*conductivity*/ )
    : _transport_mixture( new TransportMixture<CoeffType>(mixture.chem_mixture()) ),
      _wilke_mixture( new WilkeTransportMixture<CoeffType>(*_transport_mixture) ),
      _diffusion( new MixtureDiffusion<ConstantLewisDiffusivity<CoeffType>,CoeffType>(*_transport_mixture) ),
      _conductivity( new MixtureConductivity<ThermalConductivity,CoeffType>(*_transport_mixture) ),
      _micro_thermo( new typename ThermalConductivity::micro_thermo_type(mixture.chem_mixture()) ),
      _wilke_eval(new Evaluator(*_wilke_mixture,*_diffusion,viscosity,*_conductivity))
  {
    antioch_deprecated();
    Antioch::build_eucken_thermal_conductivity<typename ThermalConductivity::micro_thermo_type,CoeffType>(*_conductivity, *_micro_thermo );
  }

  template<class Viscosity, class ThermalConductivity, class CoeffType>
  WilkeEvaluator<Viscosity,ThermalConductivity,CoeffType>::~WilkeEvaluator()
  {
    delete _wilke_eval;
    delete _micro_thermo;
    delete _conductivity;
    delete _diffusion;
    delete _wilke_mixture;
    delete _transport_mixture;

    return;
  }

  template<class Viscosity, class ThermalConductivity, class CoeffType>
  template <typename VectorStateType>
  typename
    Antioch::value_type<VectorStateType>::type
        WilkeEvaluator<Viscosity,ThermalConductivity,CoeffType>::compute_phi( const VectorStateType& mu,
                     const VectorStateType& chi,
                     const unsigned int s ) const
  {

    typename Antioch::rebind<VectorStateType,VectorStateType>::type mu_mu_sqrt(mu.size());
    Antioch::init_constant(mu_mu_sqrt,mu);

    _wilke_eval->compute_mu_mu_sqrt( mu, mu_mu_sqrt);

    return _wilke_eval->compute_phi( mu_mu_sqrt, chi, s );
  }

} // end namespace Antioch

#endif // ANTIOCH_WILKE_EVALUATOR_H
