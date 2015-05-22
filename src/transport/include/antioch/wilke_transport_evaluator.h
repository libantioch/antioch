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
#include "antioch/kinetics_conditions.h"
#include "antioch/wilke_mixture.h"
#include "antioch/cmath_shims.h"
#include "antioch/diffusion_traits.h"

namespace Antioch
{


  template<class Diffusion, class Viscosity, class ThermalConductivity,
           class ThermoEvaluator,
           class Mixture,  // this is a WilkeMixture, but let's stop the templates here...
           class CoeffType=double>
  class WilkeTransportEvaluator
  {
  public:

    WilkeTransportEvaluator( const Mixture& mixture,
                             const ThermoEvaluator & thermo_evaluator,
                    const Diffusion& diffusion,
                    const Viscosity& viscosity,
                    const ThermalConductivity& conductivity );

    ~WilkeTransportEvaluator(){};

    //! mixture level diffusion, array of values
    template <typename StateType, typename VectorStateType>
    void D( const StateType & rho,
            const StateType& T,
            const VectorStateType& mass_fractions,
            VectorStateType & D_vec) const;

    //! mixture level viscosity, one value
    template <typename Conditions, typename VectorStateType>
    typename value_type<VectorStateType>::type
         mu( const Conditions& conditions,
             const VectorStateType& mass_fractions ) const;

    //! mixture level thermal conduction, one value
    template <typename Conditions, typename VectorStateType>
    typename value_type<VectorStateType>::type
        k( const Conditions & conditions,
           const VectorStateType& mass_fractions ) const;

    //! mixture level thermal conduction and viscosity, one value
    template <typename Conditions, typename StateType, typename VectorStateType>
    void mu_and_k( const Conditions& conditions,
                   const VectorStateType& mass_fractions,
                   StateType& mu, StateType& k ) const;

    //! mixture level thermal conduction, viscosity (one value) and diffusion (array of values)
    template <typename Conditions, typename StateType, typename VectorStateType>
    void mu_and_k_and_D( const Conditions& conditions, const StateType & rho,
                         const VectorStateType& mass_fractions,
                         StateType& mu, StateType& k, VectorStateType & D_vec ) const;

    //! Helper function to reduce code duplication.
    /*! Populates species viscosities and the intermediate \chi variable
        needed for Wilke's mixing rule. */
    template <typename Conditions, typename VectorStateType>
    void compute_mu_chi( const Conditions& conditions,
                         const VectorStateType& mass_fractions,
                         VectorStateType& mu,
                         VectorStateType& chi ) const;

    //! Helper function to reduce code duplication.
    /*! Computes the intermediate \phi variable needed for Wilke's mixing rule.  */
    template <typename VectorStateType>
    typename value_type<VectorStateType>::type
    compute_phi( typename rebind<VectorStateType,VectorStateType>::type & mu_mu_sqrt,
                 const VectorStateType& chi,
                 const unsigned int s ) const;

    //! Helper function to reduce code duplication.
    /*! Computes the intermediate \mu,\mu matrix needed for Wilke's mixing rule.  */
    template <typename VectorStateType>
    void compute_mu_mu_sqrt( const VectorStateType & mu,
                             typename rebind<VectorStateType,VectorStateType>::type & mu_mu_sqrt) const;

  protected:

    //! Compute species diffusion coefficients
    /*! Uses Wilke mixing rule to compute species diffusion coefficients, D_vec,
        based on the given binary diffusion matrix, D_mat. */
    template <typename VectorStateType, typename MatrixStateType>
    void diffusion_mixing_rule( const ChemicalMixture<CoeffType> & mixture,
                                const VectorStateType & mass_fractions,
                                const MatrixStateType & D_mat,
                                VectorStateType & D_vec ) const;

    const Mixture& _mixture;

    const ThermoEvaluator& _thermo;

    const Diffusion& _diffusion;

    const Viscosity& _viscosity;

    const ThermalConductivity& _conductivity;

  private:

    WilkeTransportEvaluator();

  };

  template<class Diff, class V, class TC, class ThermoEvaluator, class Mixture, class CoeffType>
  WilkeTransportEvaluator<Diff,V,TC,ThermoEvaluator,Mixture,CoeffType>::WilkeTransportEvaluator( const Mixture& mixture,
                                                                                                 const ThermoEvaluator& thermo_evaluator,
                                                                                                 const Diff& diffusion,
                                                                                                 const V& viscosity,
                                                                                                 const TC& conductivity )
    : _mixture(mixture),
      _thermo(thermo_evaluator),
      _diffusion(diffusion),
      _viscosity(viscosity),
      _conductivity(conductivity)
  {
    return;
  }

  template<class Diff, class V, class TC, class ThermoEvaluator,class Mixture, class CoeffType>
  template <typename StateType, typename VectorStateType>
  void WilkeTransportEvaluator<Diff,V,TC,ThermoEvaluator,Mixture,CoeffType>::D( const StateType & rho,
                                                                                const StateType& T,
                                                                                const VectorStateType& mass_fractions,
                                                                                VectorStateType & D_vec) const
  {
#ifdef ANTIOCH_HAVE_CXX_STATIC_ASSERT
    static_assert( !DiffusionTraits<typename Diff::Type,CoeffType>::is_binary_diffusion,
                   "This function requires a binary diffusion model to compute D!" );
#endif

    typename rebind<VectorStateType,VectorStateType>::type D_mat(D_vec.size());
    init_constant(D_mat,D_vec);

     const StateType molar_density = rho / _mixture.chem_mixture().M(mass_fractions); // total molar density

    _diffusion.compute_binary_diffusion_matrix( T, molar_density, D_mat );

    this->diffusion_mixing_rule( _mixture().chem_mixture(),
                                 mass_fractions,
                                 D_mat,
                                 D_vec );
  }

  template<class Diff, class V, class TC, class ThermoEvaluator, class Mixture, class CoeffType>
  template <typename Conditions, typename VectorStateType>
  inline
  typename value_type<VectorStateType>::type
         WilkeTransportEvaluator<Diff,V,TC,ThermoEvaluator,Mixture,CoeffType>::mu( const Conditions& conditions,
                                                                                   const VectorStateType& mass_fractions ) const
  {

    typename constructor_or_reference<const KineticsConditions<typename value_type<VectorStateType>::type,VectorStateType>, const Conditions>::type  //either (KineticsConditions<> &) or (KineticsConditions<>)
                transport_conditions(conditions);

    typename value_type<VectorStateType>::type mu_mix = zero_clone(transport_conditions.T());

    VectorStateType mu  = zero_clone(mass_fractions);
    VectorStateType chi = zero_clone(mass_fractions);


    typename Antioch::rebind<VectorStateType,VectorStateType>::type mu_mu_sqrt(mu.size());
    Antioch::init_constant(mu_mu_sqrt,mu);

    this->compute_mu_chi( transport_conditions, mass_fractions, mu, chi );

    this->compute_mu_mu_sqrt( mu, mu_mu_sqrt);

    for( unsigned int s = 0; s < _mixture.chem_mixture().n_species(); s++ )
      {
        typename value_type<VectorStateType>::type phi_s = this->compute_phi( mu_mu_sqrt, chi, s );

        mu_mix += mu[s]*chi[s]/phi_s;
      }

    return mu_mix;
  }

  template<class D, class V, class TC, class ThermoEvaluator, class Mixture, class CoeffType>
  template <typename Conditions, typename VectorStateType>
  typename value_type<VectorStateType>::type
        WilkeTransportEvaluator<D,V,TC,ThermoEvaluator,Mixture,CoeffType>::k( const Conditions& conditions,
                                                                              const VectorStateType& mass_fractions ) const
  {
    antioch_assert_equal_to(mass_fractions.size(), _mixture.chem_mixture().n_species());

    typename constructor_or_reference<const KineticsConditions<typename value_type<VectorStateType>::type,VectorStateType>, const Conditions>::type  //either (KineticsConditions<> &) or (KineticsConditions<>)
                transport_conditions(conditions);

    typename value_type<VectorStateType>::type k_mix = zero_clone(transport_conditions.T());

    VectorStateType mu  = zero_clone(mass_fractions);
    VectorStateType chi = zero_clone(mass_fractions);

    typename Antioch::rebind<VectorStateType,VectorStateType>::type mu_mu_sqrt(mu.size());
    Antioch::init_constant(mu_mu_sqrt,mu);

    this->compute_mu_chi( transport_conditions, mass_fractions, mu, chi );

    this->compute_mu_mu_sqrt( mu, mu_mu_sqrt);

    for( unsigned int s = 0; s < _mixture.chem_mixture().n_species(); s++ )
      {
        typename value_type<VectorStateType>::type  phi_s = this->compute_phi( mu_mu_sqrt, chi, s );

        typename value_type<VectorStateType>::type  k_s = zero_clone(transport_conditions.T());

        k_s =  _conductivity.conductivity_without_diffusion( s,
                                                             transport_conditions.T(),
                                                             mu[s] );

        k_mix += k_s*chi[s]/phi_s;
      }

    return k_mix;
  }

  template<class D, class V, class TC, class ThermoEvaluator, class Mixture, class CoeffType>
  template <typename Conditions, typename StateType, typename VectorStateType>
  void WilkeTransportEvaluator<D,V,TC,ThermoEvaluator,Mixture,CoeffType>::mu_and_k( const Conditions& conditions,
                                                                                    const VectorStateType& mass_fractions,
                                                                                    StateType& mu_mix,
                                                                                    StateType& k_mix ) const
  {
    typename constructor_or_reference<const KineticsConditions<StateType,VectorStateType>, const Conditions>::type  //either (KineticsConditions<> &) or (KineticsConditions<>)
                transport_conditions(conditions);

    mu_mix = zero_clone(transport_conditions.T());
    k_mix  = zero_clone(transport_conditions.T());

    VectorStateType mu  = zero_clone(mass_fractions);
    VectorStateType chi = zero_clone(mass_fractions);

    typename Antioch::rebind<VectorStateType,VectorStateType>::type mu_mu_sqrt(mu.size());
    Antioch::init_constant(mu_mu_sqrt,mu);

    this->compute_mu_chi( transport_conditions, mass_fractions, mu, chi );

    this->compute_mu_mu_sqrt( mu, mu_mu_sqrt);

    for( unsigned int s = 0; s < _mixture.transport_mixture().n_species(); s++ )
      {
        StateType phi_s = this->compute_phi( mu_mu_sqrt, chi, s );

        StateType k_s = zero_clone(transport_conditions.T());

        k_s =  _conductivity.conductivity_without_diffusion( s,
                                                             transport_conditions.T(),
                                                             mu[s] );

        mu_mix += mu[s]*chi[s]/phi_s;
        k_mix += k_s*chi[s]/phi_s;
      }

    return;
  }

  template<class Diff, class V, class TC, class ThermoEvaluator, class Mixture, class CoeffType>
  template <typename Conditions, typename StateType, typename VectorStateType>
  void WilkeTransportEvaluator<Diff,V,TC,ThermoEvaluator,Mixture,CoeffType>::mu_and_k_and_D( const Conditions& conditions,
                                                                                             const StateType & rho,
                                                                                             const VectorStateType& mass_fractions,
                                                                                             StateType& mu_mix,
                                                                                             StateType& k_mix,
                                                                                             VectorStateType & D_vec ) const
  {

    typename constructor_or_reference<const KineticsConditions<StateType,VectorStateType>, const Conditions>::type  //either (KineticsConditions<> &) or (KineticsConditions<>)
                transport_conditions(conditions);

    mu_mix = zero_clone(transport_conditions.T());
    k_mix  = zero_clone(transport_conditions.T());
    D_vec = zero_clone(mass_fractions);

    VectorStateType mu  = zero_clone(mass_fractions);
    VectorStateType k  = zero_clone(mass_fractions);
    VectorStateType chi = zero_clone(mass_fractions);

    typedef typename Antioch::rebind<VectorStateType,VectorStateType>::type MatrixStateType;

    MatrixStateType mu_mu_sqrt(mu.size());
    Antioch::init_constant(mu_mu_sqrt,mu);

    MatrixStateType D_mat(mass_fractions.size());
    init_constant(D_mat,D_vec);


    this->compute_mu_chi( transport_conditions, mass_fractions, mu, chi );

    this->compute_mu_mu_sqrt( mu, mu_mu_sqrt);

    const StateType molar_density = rho / _mixture.chem_mixture().M(mass_fractions); // total molar density

    // If we're using a binary diffusion model, compute D_mat, D_vec now
    if( DiffusionTraits<typename Diff::Type,CoeffType>::is_binary_diffusion )
      {
        _diffusion.compute_binary_diffusion_matrix(transport_conditions.T(), molar_density, D_mat);

        this->diffusion_mixing_rule<VectorStateType,MatrixStateType>( _mixture.chem_mixture(),
                                                                      mass_fractions,
                                                                      D_mat,
                                                                      D_vec );
      }

    for( unsigned int s = 0; s < _mixture.transport_mixture().n_species(); s++ )
      {
        StateType phi_s = this->compute_phi( mu_mu_sqrt, chi, s );

        k[s] = _conductivity.conductivity_with_diffusion( s,
                                                          transport_conditions.T(),
                                                          molar_density*_mixture.chem_mixture().M(s),
                                                          //rho*mass_fractions[s], // species density, rho_s
                                                          mu[s],
                                                          D_mat[s][s] );

        mu_mix += mu[s]*chi[s]/phi_s;
        k_mix += k[s]*chi[s]/phi_s;
      }

    if( DiffusionTraits<typename Diff::Type,CoeffType>::is_species_diffusion )
      {
        for( unsigned int s = 0; s < _mixture.transport_mixture().n_species(); s++ )
          {
            _diffusion.compute_species_diffusivity(s,
                                                   rho,
                                                   _thermo.cp(transport_conditions.temp_cache(),s),
                                                   k[s],
                                                   D_vec[s]);
          }
      }

    return;
  }

  template<class D, class V, class TC, class ThermoEvaluator, class Mixture, class CoeffType>
  template <typename Conditions, typename VectorStateType>
  void WilkeTransportEvaluator<D,V,TC,ThermoEvaluator,Mixture,CoeffType>::compute_mu_chi( const Conditions& conditions,
                                                                                          const VectorStateType& mass_fractions,
                                                                                          VectorStateType& mu,
                                                                                          VectorStateType& chi ) const
  {

    typename constructor_or_reference<const KineticsConditions<typename value_type<VectorStateType>::type,VectorStateType>, const Conditions>::type  //either (KineticsConditions<> &) or (KineticsConditions<>)
                transport_conditions(conditions);

    const typename value_type<VectorStateType>::type M = _mixture.chem_mixture().M(mass_fractions);

    // Precompute needed quantities
    // chi_s = w_s*M/M_s
    for( unsigned int s = 0; s < _mixture.chem_mixture().n_species(); s++ )
      {
        mu[s] = _viscosity(s,transport_conditions.T());
        chi[s] = mass_fractions[s]*M/_mixture.chem_mixture().M(s);
      }

    return;
  }

  template<class D, class V, class TC, class ThermoEvaluator, class Mixture, class CoeffType>
  template <typename VectorStateType>
  typename
  value_type<VectorStateType>::type
  WilkeTransportEvaluator<D,V,TC,ThermoEvaluator,Mixture,CoeffType>::compute_phi(typename Antioch::rebind<VectorStateType,VectorStateType>::type & mu_mu_sqrt,
                                                                                 const VectorStateType& chi,
                                                                                 const unsigned int s ) const
  {
    using std::sqrt;

    typedef typename value_type<VectorStateType>::type StateType;

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

  template<class D, class V, class TC, class ThermoEvaluator, class Mixture, class CoeffType>
  template <typename VectorStateType>
  inline
  void WilkeTransportEvaluator<D,V,TC,ThermoEvaluator,Mixture,CoeffType>::compute_mu_mu_sqrt( const VectorStateType & mu,
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

  template<class D, class V, class TC, class ThermoEvaluator, class Mixture, class CoeffType>
  template <typename VectorStateType, typename MatrixStateType>
  void WilkeTransportEvaluator<D,V,TC,ThermoEvaluator,Mixture,CoeffType>::diffusion_mixing_rule( const ChemicalMixture<CoeffType> & mixture,
                                                                                                 const VectorStateType & mass_fractions,
                                                                                                 const MatrixStateType & D_mat,
                                                                                                 VectorStateType & D_vec ) const
  {
    antioch_assert_equal_to(D_vec.size(),mixture.n_species());
    antioch_assert_equal_to(D_vec.size(),mass_fractions.size());
    antioch_assert_equal_to(D_vec.size(),D_mat.size());

#ifndef NDEBUG
    // D_Mat should be square
    for(unsigned int s = 0; s < D_vec.size(); s++)
      antioch_assert_equal_to(D_vec.size(),D_mat[s].size());
#endif

    VectorStateType molar_fractions = zero_clone(mass_fractions);

    mixture.X(mixture.M(mass_fractions),mass_fractions,molar_fractions);

    // D_s = (1 - Y_s) / (sum_{j \neq s} x_j/D_{s,j})
    for(unsigned int s = 0; s < D_vec.size(); s++)
      {
        D_vec[s] = constant_clone(mass_fractions[s],1) - mass_fractions[s];

        typename value_type<VectorStateType>::type denom = zero_clone(mass_fractions[0]);

        for(unsigned int j = 0; j < D_vec.size(); j++)
          {
            if(j == s)
              continue;

            denom += molar_fractions[j] / D_mat[s][j];
          }

        D_vec[s] /= denom;
      }
  }

} // end namespace Antioch

#endif // ANTIOCH_WILKE_TRANSPORT_EVALUATOR_H
