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

#ifndef ANTIOCH_MIXTURE_CONDUCTIVITY_H
#define ANTIOCH_MIXTURE_CONDUCTIVITY_H

#include "antioch/mixture_transport_base.h"

namespace Antioch
{
  //! Container class for species thermal conductivities
  /*! For the given set of chemical species in the input TransportMixture, this contains
    all the thermal conductivities for each of those species and provides and interface for
    computing the species thermal conductivity. Total conductivity is computed by a mixing model,
    e.g. WilkeTransportMixture. This class is templated on the conductivity model,
    so an inherent assumption is that all species conductivities have the same model. */
  template<typename Conductivity, typename MicroThermo, class CoeffType=double>
  class MixtureConductivity : public MixtureTransportBase<CoeffType>
  {
  public:

    MixtureConductivity( const TransportMixture<CoeffType>& transport_mixture );

    virtual ~MixtureConductivity();

    //! Add species viscosity
    void add( unsigned int s,
              const std::vector<CoeffType>& coeffs,
              const MicroThermo& thermo );

    //! Reset model coefficients for viscosity model of species s
    void reset_coeffs( const unsigned int s,
                       const std::vector<CoeffType> coeffs );

    template<typename StateType>
    StateType conductivity_with_diffusion( unsigned int s,
                                           const StateType& T,
                                           const StateType& rho,
                                           const StateType& mu_s,
                                           const StateType& D_ss ) const;

    template <typename StateType>
    StateType conductivity_without_diffusion( unsigned int s,
                                              const StateType& T,
                                              const StateType& mu_s ) const;

  protected:

    std::vector<Conductivity*> _species_conductivities;

  private:

    MixtureConductivity();

  };

  template<typename Conductivity, typename MicroThermo, class CoeffType>
  MixtureConductivity<Conductivity,MicroThermo,CoeffType>::MixtureConductivity(  const TransportMixture<CoeffType>& transport_mixture )
    :  MixtureTransportBase<CoeffType>(transport_mixture),
    _species_conductivities( transport_mixture.n_species(), NULL )
  {}

  template<typename Conductivity, typename MicroThermo, class CoeffType>
  MixtureConductivity<Conductivity,MicroThermo,CoeffType>::~MixtureConductivity()
  {
    // Need to delete all the species viscosities we allocated
    for( typename std::vector<Conductivity*>::iterator it = _species_conductivities.begin();
         it != _species_conductivities.end(); ++it )
      {
        delete (*it);
      }
  }

  template<typename Conductivity, typename MicroThermo, class CoeffType>
  void MixtureConductivity<Conductivity,MicroThermo,CoeffType>::add( unsigned int s,
                                                                     const std::vector<CoeffType>& coeffs,
                                                                     const MicroThermo& thermo )
  {
    antioch_assert_less_equal( s, _species_conductivities.size() );
    antioch_assert( !_species_conductivities[s] );

    _species_conductivities[s] = new Conductivity(thermo, coeffs);
  }

  template<typename Conductivity, typename MicroThermo, class CoeffType>
  void MixtureConductivity<Conductivity,MicroThermo,CoeffType>::reset_coeffs( const unsigned int s,
                                                                              const std::vector<CoeffType> coeffs )
  {
    _species_conductivities[s]->reset_coeffs(coeffs);
  }

  template<typename Conductivity, typename MicroThermo, class CoeffType>
  template<typename StateType>
  StateType MixtureConductivity<Conductivity,MicroThermo,CoeffType>::conductivity_with_diffusion( unsigned int s,
                                                                                                  const StateType& T,
                                                                                                  const StateType& rho,
                                                                                                  const StateType& mu_s,
                                                                                                  const StateType& D_ss ) const
  {
    return (*this->_species_conductivities[s])(s,mu_s,T,rho,D_ss);
  }

  template<typename Conductivity, typename MicroThermo, class CoeffType>
  template<typename StateType>
  StateType MixtureConductivity<Conductivity,MicroThermo,CoeffType>::conductivity_without_diffusion( unsigned int s,
                                                                                                     const StateType& T,
                                                                                                     const StateType& mu_s ) const
  {
    return (*this->_species_conductivities[s])(s,mu_s,T);
  }

} // end namespace Antioch

#endif // ANTIOCH_MIXTURE_CONDUCTIVITY_H
