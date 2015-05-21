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

#ifndef ANTIOCH_MIXTURE_DIFFUSION_BASE_H
#define ANTIOCH_MIXTURE_DIFFUSION_BASE_H

#include "antioch/mixture_transport_base.h"

namespace Antioch
{

  template<typename Subclass, class CoeffType>
  class MixtureDiffusionBase : public MixtureTransportBase<CoeffType>
  {
  public:

    MixtureDiffusionBase( const TransportMixture<CoeffType>& transport_mixture );

    virtual ~MixtureDiffusionBase(){};

    template<typename StateType, typename MatrixStateType>
    void compute_binary_diffusion_matrix( const StateType& T, const StateType& molar_density, MatrixStateType& D ) const;

    //! Computes species diffusivity for species s
    template<typename StateType>
    void compute_species_diffusivity( unsigned int s, const StateType& rho,
                                      const StateType& cp, const StateType& k,
                                      StateType& D ) const;
  };

  template<typename Subclass, class CoeffType>
  MixtureDiffusionBase<Subclass,CoeffType>::MixtureDiffusionBase( const TransportMixture<CoeffType>& transport_mixture )
    :  MixtureTransportBase<CoeffType>(transport_mixture)
  {}

  template<typename Subclass, class CoeffType>
  template<typename StateType, typename MatrixStateType>
  void MixtureDiffusionBase<Subclass,CoeffType>::compute_binary_diffusion_matrix( const StateType& T,
                                                                                  const StateType& molar_density,
                                                                                  MatrixStateType& D ) const
  {
    static_cast<const Subclass*>(this)->compute_binary_diffusion_matrix_impl(T,molar_density,D);
  }

  template<typename Subclass, class CoeffType>
  template<typename StateType>
  void MixtureDiffusionBase<Subclass,CoeffType>::compute_species_diffusivity( unsigned int s,
                                                                              const StateType& rho,
                                                                              const StateType& cp,
                                                                              const StateType& k,
                                                                              StateType& D ) const
  {
    static_cast<const Subclass*>(this)->compute_species_diffusivity_impl(s,rho,cp,k,D);
  }

} // end namespace Antioch

#endif // ANTIOCH_MIXTURE_DIFFUSION_BASE_H
