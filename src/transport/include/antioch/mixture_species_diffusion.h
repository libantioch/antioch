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

#ifndef ANTIOCH_MIXTURE_SPECIES_DIFFUSION_H
#define ANTIOCH_MIXTURE_SPECIES_DIFFUSION_H

#include "antioch/mixture_transport_base.h"

namespace Antioch
{
  //! Container class for species binary diffusion models
  /*! For the given set of chemical species in the input TransportMixture, this contains
   *  the diffusion models for each of those species and provides and interface for
   *  computing the species diffusion coefficients. Total diffusion coefficients is computed
   *  by a mixing model, e.g. WilkeTransportMixture. This class is templated on the diffusion model,
   *  so an inherent assumption is that all species diffusions have the same model.
   */
  template<typename Diffusion, class CoeffType=double>
  class MixtureSpeciesDiffusion : public MixtureDiffusionBase<MixtureSpeciesDiffusion<Diffusion,CoeffType>,CoeffType>
  {
  public:

    MixtureSpeciesDiffusion( const TransportMixture<CoeffType>& transport_mixture );

    virtual ~MixtureSpeciesDiffusion();

    //! Add species diffusion model for species s
    void add( unsigned int s, const std::vector<CoeffType>& coeffs );

    //! Reset model coefficients for species diffusion model for species s
    void reset_coeffs( const unsigned int s, const std::vector<CoeffType> coeffs );

    //! Define the diffusion model. Mainly for use in DiffusionTraits
    typedef Diffusion Type;

    //! Friend the base class so we can make the implementation protected
    friend class MixtureDiffusionBase<MixtureBinaryDiffusion<Diffusion,CoeffType>,CoeffType>;

  protected:

    //! Stores diffusivity models
    std::vector<Diffusion*> _species_diffusivities;

    //! Computes species diffusivity for species s
    template<typename StateType>
    void compute_species_diffusivity_impl( unsigned int s, const StateType& rho,
                                           const StateType& cp, const StateType& k,
                                           StateType& D ) const;

    //! Should never be called by this class, facilitating interaction in WilkeTransportEvaluator
    /*!
     * There are wildly differing requirements between species and binary diffusion
     * models, but we want to have a uniform interface through WilkeTransportEvaluator.
     * So we need to implement this function, but it will throw an error.
     */
    template<typename StateType, typename MatrixStateType>
    void compute_binary_diffusion_matrix_impl( const StateType& T,
                                               const StateType& molar_density,
                                               MatrixStateType& D ) const;

  private:

    MixtureSpeciesDiffusion();
  };

  template<typename Diffusion, class CoeffType>
  MixtureSpeciesDiffusion<Diffusion,CoeffType>::MixtureSpeciesDiffusion( const TransportMixture<CoeffType>& transport_mixture )
    :  MixtureDiffusionBase<MixtureSpeciesDiffusion<Diffusion,CoeffType>,CoeffType>(transport_mixture),
    _species_diffusivities( transport_mixture.n_species(), NULL )
  {
    //#ifdef ANTIOCH_HAVE_CXX_STATIC_ASSERT
    static_assert( DiffusionTraits<Diffusion,CoeffType>::is_species_diffusion,
                   "Can only instantiate MixtureSpeciesDiffusion with a species diffusion model!" );
    //#endif
  }

  template<typename Diffusion, class CoeffType>
  MixtureSpeciesDiffusion<Diffusion,CoeffType>::~MixtureSpeciesDiffusion()
  {
    for( typename std::vector<Diffusion*>::iterator it = _species_diffusivities.begin();
         it != _species_diffusivities.end(); ++it )
      delete *it;
  }

  template<typename Diffusion, class CoeffType>
  void MixtureSpeciesDiffusion<Diffusion,CoeffType>::add( unsigned int s,
                                                          const std::vector<CoeffType>& coeffs )
  {
    antioch_assert_less( s, _species_diffusivities.size() );

    _species_diffusivities[s] = new Diffusion( coeffs );
  }

  template<typename Diffusion, class CoeffType>
  void MixtureSpeciesDiffusion<Diffusion,CoeffType>::reset_coeffs( const unsigned int s,
                                                                   const std::vector<CoeffType> coeffs )
  {
    antioch_assert_less( s, _species_diffusivities.size() );
    antioch_assert(_species_diffusivities[s]);

    _species_diffusivities[s]->reset_coeffs( coeffs );
  }

  template<typename Diffusion, class CoeffType>
  template<typename StateType>
  inline
  void MixtureSpeciesDiffusion<Diffusion,CoeffType>::compute_species_diffusivity_impl( unsigned int s,
                                                                                       const StateType& rho,
                                                                                       const StateType& cp,
                                                                                       const StateType& k,
                                                                                       StateType& D ) const
  {
    antioch_assert_less( s, _species_diffusivities.size());
    antioch_assert(_species_diffusivities[s]);

    (*_species_diffusivities[s])(rho,cp,k,D);
  }

  template<typename Diffusion, class CoeffType>
  template<typename StateType, typename MatrixStateType>
  inline
  void MixtureSpeciesDiffusion<Diffusion,CoeffType>::compute_binary_diffusion_matrix_impl( const StateType& T,
                                                                                           const StateType& molar_density,
                                                                                           MatrixStateType& D ) const
  {
    std::string error = "ERROR: You're trying to use a binary diffusion implementation\n";
    error += "       with a species diffusion model!\n";
    antioch_msg_error(error);
  }

} // end namespace Antioch

#endif // ANTIOCH_MIXTURE_SPECIES_DIFFUSION_H
