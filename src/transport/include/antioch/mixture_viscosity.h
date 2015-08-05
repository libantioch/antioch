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

#ifndef ANTIOCH_MIXTURE_VISCOSITY_H
#define ANTIOCH_MIXTURE_VISCOSITY_H

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/mixture_transport_base.h"
#include "antioch/species_viscosity_base.h"
// C++
#include <string>
#include <vector>

namespace Antioch
{
  //! Container class for species viscosities
  /*! For the given set of chemical species in the input TransportMixture, this contains
      all the viscosities for each of those species and provides and interface for
      computing the species viscosity. Total viscosity is computed by a mixing model,
      e.g. MixtureAveragedTransportMixture. This class is templated on the viscosity model,
      so an inherent assumption is that all species viscosities have the same model. */
  template<typename Viscosity, class CoeffType=double>
  class MixtureViscosity : public MixtureTransportBase<CoeffType>
  {
  public:

    MixtureViscosity( const TransportMixture<CoeffType>& transport_mixture );
    ~MixtureViscosity();

    //! Evaluate viscosity for species s
    /*! Total viscosity computed by mixing model, e.g. MixtureAveragedTransportEvaluator */
    template <typename StateType>
    StateType operator()( const unsigned int s, const StateType& T ) const;

    //! Add species viscosity
    void add( const std::string& species_name,
	      const std::vector<CoeffType>& coeffs );

    //! Reset model coefficients for viscosity model of species s
    void reset_coeffs( const unsigned int s,
                       const std::vector<CoeffType> coeffs );

    const std::vector<Viscosity*> & species_viscosities() const;

    //! Formatted print, by default to \p std::cout
    void print(std::ostream& os = std::cout) const;

    //! Formatted print.
    friend std::ostream& operator<<(std::ostream& os, const MixtureViscosity& mu)
    {
      mu.print(os);
      return os;
    }

    typedef Viscosity species_viscosity_type;

  protected:

    std::vector<SpeciesViscosityBase<Viscosity,CoeffType>*> _species_viscosities;

  private:

    MixtureViscosity();

  };

  template<typename Viscosity, class CoeffType>
  MixtureViscosity<Viscosity,CoeffType>::MixtureViscosity(  const TransportMixture<CoeffType>& transport_mixture )
    :  MixtureTransportBase<CoeffType>(transport_mixture),
       _species_viscosities( transport_mixture.n_species(), NULL )
  {}

  template<typename Viscosity, class CoeffType>
  MixtureViscosity<Viscosity,CoeffType>::~MixtureViscosity()
  {
    // Need to delete all the species viscosities we allocated
    for( typename std::vector<SpeciesViscosityBase<Viscosity,CoeffType>*>::iterator it = _species_viscosities.begin();
	 it != _species_viscosities.end(); ++it )
      {
	delete (*it);
      }
  }

  template<typename Viscosity, class CoeffType>
  void MixtureViscosity<Viscosity,CoeffType>::add( const std::string& species_name,
						   const std::vector<CoeffType>& coeffs )
  {
    antioch_assert( this->_transport_mixture.species_name_map().find(species_name) !=
		    this->_transport_mixture.species_name_map().end() );

    unsigned int s = this->_transport_mixture.species_name_map().find(species_name)->second;

    antioch_assert_less_equal( s, _species_viscosities.size() );
    antioch_assert( !_species_viscosities[s] );

    _species_viscosities[s] = new Viscosity(coeffs);
    return;
  }

  template<typename Viscosity, class CoeffType>
  void MixtureViscosity<Viscosity,CoeffType>::reset_coeffs( const unsigned int s,
                                                            const std::vector<CoeffType> coeffs )
  {
    _species_viscosities[s]->reset_coeffs(coeffs);
  }

  template<typename Viscosity, class CoeffType>
  template<typename StateType>
  inline
  StateType MixtureViscosity<Viscosity,CoeffType>::operator()( const unsigned int s,
							       const StateType& T ) const
  {
    antioch_assert_less_equal( s, _species_viscosities.size() );
    antioch_assert( _species_viscosities[s] );

    return (*_species_viscosities[s])(T);
  }

  template<typename Viscosity, class CoeffType>
  inline
  const std::vector<Viscosity*>& MixtureViscosity<Viscosity,CoeffType>::species_viscosities() const
  {
    return _species_viscosities;
  }

  template<typename Viscosity, class CoeffType>
  void MixtureViscosity<Viscosity,CoeffType>::print(std::ostream& os) const
  {
    antioch_assert_equal_to( _species_viscosities.size(), this->_transport_mixture.n_species() );

    for( unsigned int s = 0; s < this->_transport_mixture.n_species(); s++ )
      {
	const Species& species = this->_transport_mixture.species_list()[s];
	const std::string& name = this->_transport_mixture.species_inverse_name_map().find( species  )->second;

	os << "mu(" << name << ") = " << (*_species_viscosities[s]) << std::endl;
      }

    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_MIXTURE_VISCOSITY_H
