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
#include "antioch/mixture_averaged_transport_mixture.h"
#include "antioch/cmath_shims.h"
#include "antioch/mixture_diffusion.h"
#include "antioch/mixture_viscosity.h"
#include "antioch/mixture_conductivity.h"
#include "antioch/diffusion_traits.h"
#include "antioch/conductivity_traits.h"

namespace Antioch
{
  //! Compute transport properties using ``mixture averaged" model
  /*!
   * Use the species transport models (template parameters) to evaluate species
   * values, then use Wilke's mixing rule to compute the mixture viscosity and thermal conductivity
   * and the appropriate DiffusivityType for the mixture diffusivity. This is
   * the expected interface for the user. Underlying compile time decisions
   * are made based on the species transport models that should be invisible
   * to the user. The preferred, most efficient, and most general
   * method is mu_and_k_and_D for evaluating the transport quantities.
   * The other methods are less efficient and may only be available
   * for a subset of species model; they are present for backwards
   * compatibility and pedagogical reasons.
   *
   * In a multi-threaded computing environment, this object should be created *after* threads
   * have been forked. In other words, this object is meant to live only temporarily for computing
   * quantities.
   */
  template<class Diffusion, class Viscosity, class ThermalConductivity,class CoeffType=double>
  class MixtureAveragedTransportEvaluator
  {
  public:

    MixtureAveragedTransportEvaluator( const MixtureAveragedTransportMixture<CoeffType>& mixture,
                                       const MixtureDiffusion<Diffusion,CoeffType>& diffusion,
                                       const MixtureViscosity<Viscosity,CoeffType>& viscosity,
                                       const MixtureConductivity<ThermalConductivity,CoeffType>& conductivity );

    ~MixtureAveragedTransportEvaluator(){};

    /*!
     *  For the mixture averaged diffusion models, there are three typical ways
     *  of computing an average diffusivity, depending on the "units" of the flux
     *  and which variable with which the gradient is taken:
     *  -# molar flux with gradients of mole fraction (MOLE_FLUX_MOLE_FRACTION)
     *  -# mass flux with gradients of mass fraction (MASS_FLUX_MASS_FRACTION)
     *  -# mass flux with gradients of mole fraction (MASS_FLUX_MOLE_FRACTION)
     *  See, for example: Kee, Coltrin, Glarborg, ``Chemically Reacting Flow", Wiley, 2003, Chapter 12.7
     */
    enum DiffusivityType { MOLE_FLUX_MOLE_FRACTION,
                           MASS_FLUX_MASS_FRACTION,
                           MASS_FLUX_MOLE_FRACTION };

    //! Mixture diffusivities for each species, in [m^2/s]
    /*! Only valid for BinaryDiffusionBase species diffusion models.
     *  Compile time error if otherwise.
     */
    template <typename StateType, typename VectorStateType>
    void D( const StateType & rho,
            const StateType& T,
            const VectorStateType& mass_fractions,
            VectorStateType & D_vec,
            DiffusivityType diff_type = DiffusivityType::MASS_FLUX_MOLE_FRACTION ) const;

    //! Mixture viscosity, in [Pa-s]
    template <typename StateType, typename VectorStateType>
    typename value_type<VectorStateType>::type
    mu( const StateType& T, const VectorStateType& mass_fractions ) const;

    //! Mixture conducivity, in [W/m-K]
    /*! Only valid for "no diffusion" conductivity models.
     *  Compile time error if otherwise. */
    template <typename StateType, typename VectorStateType>
    typename value_type<VectorStateType>::type
    k( const StateType & T, const VectorStateType& mass_fractions ) const;

    //! Mixture viscosity and thermal conductivity, in [Pa-s], [W/m-K] respectively
    /*! Only valid for "no diffusion" conductivity models.
     *  Compile time error if otherwise. */
    template <typename StateType, typename VectorStateType>
    void mu_and_k( const StateType& T, const VectorStateType& mass_fractions,
                   StateType& mu, StateType& k ) const;

    //! Mixture viscosity, thermal conductivity, and diffusivities in [Pa-s], [W/m-K], [m^2/s] respectively
    /*! This is the preferred, most efficient, and most general method. */
    template <typename StateType, typename VectorStateType>
    void mu_and_k_and_D( const StateType& T, const StateType& rho, const StateType& cp,
                         const VectorStateType& mass_fractions,
                         StateType& mu, StateType& k, VectorStateType& D_vec,
                         DiffusivityType diff_type = DiffusivityType::MASS_FLUX_MOLE_FRACTION ) const;

    //! Helper function to reduce code duplication.
    /*! Populates species viscosities and the intermediate \chi variable
     *  needed for Wilke's mixing rule. This is not intended for the
     *  user; only public to facilitate testing.
     */
    template <typename StateType, typename VectorStateType>
    void compute_mu_chi( const StateType& T,
                         const VectorStateType& mass_fractions,
                         VectorStateType& mu,
                         VectorStateType& chi ) const;

    //! Helper function to reduce code duplication.
    /*! Computes the intermediate \phi variable needed for Wilke's mixing rule.
     *  Not intended for the user; only public to facilitate testing.
     */
    template <typename VectorStateType>
    typename value_type<VectorStateType>::type
    compute_phi( typename rebind<VectorStateType,VectorStateType>::type & mu_mu_sqrt,
                 const VectorStateType& chi,
                 const unsigned int s ) const;

    //! Helper function to reduce code duplication.
    /*! Computes the intermediate \mu,\mu matrix needed for Wilke's mixing rule.
     *  Not intended for the user; only public to facilitate testing.
     */
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
                                DiffusivityType diff_type,
                                VectorStateType & D_vec ) const;

    const MixtureAveragedTransportMixture<CoeffType>& _mixture;

    const MixtureDiffusion<Diffusion,CoeffType>& _diffusion;

    const MixtureViscosity<Viscosity,CoeffType>& _viscosity;

    const MixtureConductivity<ThermalConductivity,CoeffType>& _conductivity;

  private:

    MixtureAveragedTransportEvaluator();

  };

  template<class Diff, class Visc, class TherCond, class CoeffType>
  MixtureAveragedTransportEvaluator<Diff,Visc,TherCond,CoeffType>::MixtureAveragedTransportEvaluator( const MixtureAveragedTransportMixture<CoeffType>& mixture,
                                                                                                      const MixtureDiffusion<Diff,CoeffType>& diffusion,
                                                                                                      const MixtureViscosity<Visc,CoeffType>& viscosity,
                                                                                                      const MixtureConductivity<TherCond,CoeffType>& conductivity )
  : _mixture(mixture),
    _diffusion(diffusion),
    _viscosity(viscosity),
    _conductivity(conductivity)
  {
    return;
  }

  template<class Diff, class Visc, class TherCond, class CoeffType>
  template <typename StateType, typename VectorStateType>
  void MixtureAveragedTransportEvaluator<Diff,Visc,TherCond,CoeffType>::D( const StateType & rho,
                                                                           const StateType& T,
                                                                           const VectorStateType& mass_fractions,
                                                                           VectorStateType & D_vec,
                                                                           DiffusivityType diff_type ) const
  {
    antioch_static_assert_runtime_fallback( DiffusionTraits<Diff>::is_binary_diffusion,
                                            "ERROR: This function requires a binary diffusion model to compute D!");

    typename rebind<VectorStateType,VectorStateType>::type D_mat(D_vec.size());
    init_constant(D_mat,D_vec);

    const StateType molar_density = rho / _mixture.chem_mixture().M(mass_fractions); // total molar density

    _diffusion.compute_binary_diffusion_matrix( T, molar_density, D_mat );

    this->diffusion_mixing_rule( _mixture.chem_mixture(),
                                 mass_fractions,
                                 D_mat,
                                 diff_type,
                                 D_vec );
  }

  template<class Diff, class Visc, class TherCond, class CoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  typename value_type<VectorStateType>::type
  MixtureAveragedTransportEvaluator<Diff,Visc,TherCond,CoeffType>::mu( const StateType& T,
                                                                       const VectorStateType& mass_fractions ) const
  {
    typename value_type<VectorStateType>::type mu_mix = zero_clone(T);

    VectorStateType mu  = zero_clone(mass_fractions);
    VectorStateType chi = zero_clone(mass_fractions);


    typename Antioch::rebind<VectorStateType,VectorStateType>::type mu_mu_sqrt(mu.size());
    Antioch::init_constant(mu_mu_sqrt,mu);

    this->compute_mu_chi( T, mass_fractions, mu, chi );

    this->compute_mu_mu_sqrt( mu, mu_mu_sqrt);

    for( unsigned int s = 0; s < _mixture.chem_mixture().n_species(); s++ )
      {
        typename value_type<VectorStateType>::type phi_s = this->compute_phi( mu_mu_sqrt, chi, s );

        mu_mix += mu[s]*chi[s]/phi_s;
      }

    return mu_mix;
  }

  template<class Diff, class Visc, class TherCond, class CoeffType>
  template <typename StateType, typename VectorStateType>
  typename value_type<VectorStateType>::type
  MixtureAveragedTransportEvaluator<Diff,Visc,TherCond,CoeffType>::k( const StateType& T,
                                                                      const VectorStateType& mass_fractions ) const
  {
    antioch_static_assert_runtime_fallback( !ConductivityTraits<TherCond>::requires_diffusion,
                                            "This function requires a conductivity model independent of diffusion!");

    antioch_assert_equal_to(mass_fractions.size(), _mixture.chem_mixture().n_species());

    typename value_type<VectorStateType>::type k_mix = zero_clone(T);

    VectorStateType mu  = zero_clone(mass_fractions);
    VectorStateType chi = zero_clone(mass_fractions);

    typename Antioch::rebind<VectorStateType,VectorStateType>::type mu_mu_sqrt(mu.size());
    Antioch::init_constant(mu_mu_sqrt,mu);

    this->compute_mu_chi( T, mass_fractions, mu, chi );

    this->compute_mu_mu_sqrt( mu, mu_mu_sqrt);

    for( unsigned int s = 0; s < _mixture.chem_mixture().n_species(); s++ )
      {
        typename value_type<VectorStateType>::type  phi_s = this->compute_phi( mu_mu_sqrt, chi, s );

        typename value_type<VectorStateType>::type  k_s = zero_clone(T);

        k_s =  _conductivity.conductivity_without_diffusion( s,
                                                             T,
                                                             mu[s] );

        k_mix += k_s*chi[s]/phi_s;
      }

    return k_mix;
  }

  template<class Diff, class Visc, class TherCond, class CoeffType>
  template <typename StateType, typename VectorStateType>
  void MixtureAveragedTransportEvaluator<Diff,Visc,TherCond,CoeffType>::mu_and_k( const StateType& T,
                                                                                  const VectorStateType& mass_fractions,
                                                                                  StateType& mu_mix,
                                                                                  StateType& k_mix ) const
  {
    antioch_static_assert_runtime_fallback( !ConductivityTraits<TherCond>::requires_diffusion,
                                            "This function requires a conductivity model independent of diffusion!");

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

        k_s =  _conductivity.conductivity_without_diffusion( s, T, mu[s] );

        mu_mix += mu[s]*chi[s]/phi_s;
        k_mix += k_s*chi[s]/phi_s;
      }

    return;
  }

  template<class Diff, class Visc, class TherCond, class CoeffType>
  template <typename StateType, typename VectorStateType>
  void MixtureAveragedTransportEvaluator<Diff,Visc,TherCond,CoeffType>::mu_and_k_and_D( const StateType& T,
                                                                                        const StateType & rho,
                                                                                        const StateType& cp,
                                                                                        const VectorStateType& mass_fractions,
                                                                                        StateType& mu_mix,
                                                                                        StateType& k_mix,
                                                                                        VectorStateType & D_vec,
                                                                                        DiffusivityType diff_type ) const
  {
    antioch_static_assert_runtime_fallback( (ConductivityTraits<TherCond>::requires_diffusion &&
                                             DiffusionTraits<Diff>::is_binary_diffusion) ||
                                            (!ConductivityTraits<TherCond>::requires_diffusion &&
                                             DiffusionTraits<Diff>::is_species_diffusion),
                                            "Incompatible thermal conductivity and diffusion models!" );

    mu_mix = zero_clone(T);
    k_mix  = zero_clone(T);
    D_vec = zero_clone(mass_fractions);

    VectorStateType mu  = zero_clone(mass_fractions);
    VectorStateType k  = zero_clone(mass_fractions);
    VectorStateType chi = zero_clone(mass_fractions);

    typedef typename Antioch::rebind<VectorStateType,VectorStateType>::type MatrixStateType;

    MatrixStateType mu_mu_sqrt(mu.size());
    Antioch::init_constant(mu_mu_sqrt,mu);

    MatrixStateType D_mat(mass_fractions.size());
    init_constant(D_mat,D_vec);


    this->compute_mu_chi( T, mass_fractions, mu, chi );

    this->compute_mu_mu_sqrt( mu, mu_mu_sqrt);

    const StateType molar_density = rho / _mixture.chem_mixture().M(mass_fractions); // total molar density

    // If we're using a binary diffusion model, compute D_mat, D_vec now
    if( DiffusionTraits<Diff>::is_binary_diffusion )
      {
        _diffusion.compute_binary_diffusion_matrix(T, molar_density, D_mat);

        this->diffusion_mixing_rule<VectorStateType,MatrixStateType>( _mixture.chem_mixture(),
                                                                      mass_fractions,
                                                                      D_mat,
                                                                      diff_type,
                                                                      D_vec );
      }

    for( unsigned int s = 0; s < _mixture.transport_mixture().n_species(); s++ )
      {
        StateType phi_s = this->compute_phi( mu_mu_sqrt, chi, s );

        if( ConductivityTraits<TherCond>::requires_diffusion )
          {
            k[s] = _conductivity.conductivity_with_diffusion( s,
                                                              T,
                                                              molar_density*_mixture.chem_mixture().M(s),
                                                              //rho*mass_fractions[s], see #146
                                                              mu[s],
                                                              D_mat[s][s] );

          }
        else
          {
            k[s] =  _conductivity.conductivity_without_diffusion( s,
                                                                  T,
                                                                  mu[s] );
          }

        k_mix += k[s]*chi[s]/phi_s;
        mu_mix += mu[s]*chi[s]/phi_s;

      }



    if( DiffusionTraits<Diff>::is_species_diffusion )
      {
        for( unsigned int s = 0; s < _mixture.transport_mixture().n_species(); s++ )
          {
            _diffusion.compute_species_diffusivity(s,
                                                   rho,
                                                   cp,
                                                   k_mix,
                                                   D_vec[s]);
          }
      }

    return;
  }

  template<class Diff, class Visc, class TherCond, class CoeffType>
  template <typename StateType, typename VectorStateType>
  void MixtureAveragedTransportEvaluator<Diff,Visc,TherCond,CoeffType>::compute_mu_chi( const StateType& T,
                                                                                        const VectorStateType& mass_fractions,
                                                                                        VectorStateType& mu,
                                                                                        VectorStateType& chi ) const
  {
    const typename value_type<VectorStateType>::type M = _mixture.chem_mixture().M(mass_fractions);

    // Precompute needed quantities
    // chi_s = w_s*M/M_s
    for( unsigned int s = 0; s < _mixture.chem_mixture().n_species(); s++ )
      {
        mu[s] = _viscosity(s,T);
        chi[s] = mass_fractions[s]*M/_mixture.chem_mixture().M(s);
      }

    return;
  }

  template<class Diff, class Visc, class TherCond, class CoeffType>
  template <typename VectorStateType>
  typename
  value_type<VectorStateType>::type
  MixtureAveragedTransportEvaluator<Diff,Visc,TherCond,CoeffType>::compute_phi(typename Antioch::rebind<VectorStateType,VectorStateType>::type & mu_mu_sqrt,
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

  template<class Diff, class Visc, class TherCond, class CoeffType>
  template <typename VectorStateType>
  inline
  void MixtureAveragedTransportEvaluator<Diff,Visc,TherCond,CoeffType>::compute_mu_mu_sqrt( const VectorStateType & mu,
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

  template<class Diff, class Visc, class TherCond, class CoeffType>
  template <typename VectorStateType, typename MatrixStateType>
  void MixtureAveragedTransportEvaluator<Diff,Visc,TherCond,CoeffType>::diffusion_mixing_rule( const ChemicalMixture<CoeffType> & mixture,
                                                                                               const VectorStateType & mass_fractions,
                                                                                               const MatrixStateType & D_mat,
                                                                                               DiffusivityType diff_type,
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

    switch(diff_type)
      {
      case(MASS_FLUX_MOLE_FRACTION):
        {
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
          break;
        }
      case(MOLE_FLUX_MOLE_FRACTION):
        {
          antioch_not_implemented();
        }
      case(MASS_FLUX_MASS_FRACTION):
        {
          VectorStateType molar_fractions = zero_clone(mass_fractions);

          mixture.X(mixture.M(mass_fractions),mass_fractions,molar_fractions);

          typename value_type<VectorStateType>::type one = constant_clone(mass_fractions[0],1);

          //               term1               term2
          // 1/D_s = (sum_{j\ne s} X_j/D_{s,j}) + X_s/(1-Y_s)\sum_{j\ne s} Y_j/D_{s,j}
          for(unsigned int s = 0; s < D_vec.size(); s++)
            {
              typename value_type<VectorStateType>::type term1 = zero_clone(mass_fractions[0]);
              typename value_type<VectorStateType>::type term2 = zero_clone(mass_fractions[0]);

              for(unsigned int j = 0; j < D_vec.size(); j++)
                {
                  if(j == s)
                    continue;

                  term1 += molar_fractions[j]/D_mat[s][j];

                  term2 += mass_fractions[j]/D_mat[s][j];
                }

              term2 *=  molar_fractions[s]/(one - mass_fractions[s]);

              D_vec[s] = one/(term1+term2);
            }
          break;
        }
      default:
        {
          antioch_msg_error("ERROR: Invalid DiffusivityType in MixtureAveragedTransportEvaluator::diffusion_mixing_rule");
          antioch_error();
        }
      } // switch(diff_type)
  }

} // end namespace Antioch

#endif // ANTIOCH_WILKE_TRANSPORT_EVALUATOR_H
