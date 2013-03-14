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

#ifndef ANTIOCH_WILKE_EVALUATOR_H
#define ANTIOCH_WILKE_EVALUATOR_H

namespace Antioch
{


  template<class Viscosity, class ThermalConductivity, class CoeffType=double>
  class WilkeEvaluator
  {
  public:

    WilkeEvaluator( const WilkeMixture<CoeffType>& mixture,
		    const MixtureViscosity<CoeffType,Viscosity>& viscosity,
		    const MixtureThermalConductivity<CoeffType,ThermalConductivity>& conductivity );

    ~WilkeEvaluator();

    template <typename StateType>
    StateType mu( const StateType T,
		  const std::vector<StateType>& mass_fractions );

    template <typename StateType>
    StateType k( const StateType T,
		 const std::vector<StateType>& mass_fractions );

    template <typename StateType>
    void mu_and_k( const StateType T,
		   const std::vector<StateType>& mass_fractions,
		   StateType& mu, StateType& k );

  protected:

    const WilkeMixture<CoeffType>& _mixture;

  private:

    WilkeEvaluator();

  };

  template<class V, class T, class CoeffType>
  WilkeEvaluator<V,T,CoeffType>::WilkeEvaluator( const WilkeMixture<CoeffType>& mixture,
						 const MixtureViscosity<CoeffType,V>& viscosity,
						 const MixtureThermalConductivity<CoeffType,T>& conductivity )
    : _mixture(mixture),
      _viscosity(viscosity),
      _conductivity(conductivity)
  {
    return;
  }

  template<class V, class T, class CoeffType>
  WilkeEvaluator<CoeffType,V,T>::~WilkeEvaluator()
  {
    return;
  }

  template<class V, class T, class CoeffType>
  template <typename StateType>
  StateType WilkeEvaluator<CoeffType,V,T>::mu( const StateType T,
					       const std::vector<StateType>& mass_fractions )
  {
    StateType mu = StateType(0.0);
    
    std::vector<StateType> mu(this->n_species,0.0);
    std::vector<StateType> chi(this->n_species,0.0);
    
    const StateType M = _viscosity.chemical_mixture().M(mass_fractions);

    // Precompute needed quantities
    // chi_s = w_s*M/M_s
    for( unsigned int s = 0; s < this->n_species(); s++ )
      {
	mu[s] = _viscosity(s,T);
	chi[s] = mass_fractions[s]*M/_viscosity.chemical_mixture().M(s);
      }

    for( unsigned int s = 0; s < this->n_species(); s++ )
      {
	/* We initialize to the first iterate and loop starting from 1
	   since some StateTypes have a hard time initializing from
	   a constant. */
	// phi_s = sum_r (chi_r*(1+sqrt(mu_s/mu_r)*(Mr/Ms)^(1/4))^2)/(8*(1+Ms/Mr))
	const StateType dummy = 1.0 + std::sqrt(mu[s]/mu[0])*_mixture.Mr_Ms_to_the_one_fourth(0,s);
	StateType phi_s = chi[0]*numerator*numerator/_mixture.denom(0,s);

	for(unsigned int r = 1; r < this->n_species(); r++ )
	  {
	    const StateType numerator = 1.0 + std::sqrt(mu[s]/mu[r])*_mixture.Mr_Ms_to_the_one_fourth(r,s);
	    phi_s += chi[r]*numerator*numerator/_mixture.denom(r,s);
	  }
	
	// Now compute phi_s, chi_s
	mu += mu[s]*chi[s]/phi_s
      }

    return mu;
  }

} // end namespace Antioch

#endif // ANTIOCH_WILKE_EVALUATOR_H
