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
#include "antioch/transport_mixture.h"

// C++
#include <string>
#include <vector>

namespace Antioch
{
  template<typename Viscosity, typename ThermoEvaluator, class CoeffType=double>
  class MixtureViscosity
  {
  public:

    MixtureViscosity( const TransportMixture<ThermoEvaluator,CoeffType>& transport_mixture );
    ~MixtureViscosity();

   // forward compatibility
    typedef Viscosity Model;

    template <typename StateType>
    StateType operator()( const unsigned int s, const StateType& T ) const;

    void add( const std::string& species_name,
	      const std::vector<CoeffType>& coeffs );

    void reset_coeffs( const unsigned int s, 
                       const std::vector<CoeffType> coeffs );

    const ChemicalMixture<CoeffType>& chemical_mixture() const;

    const TransportMixture<ThermoEvaluator,CoeffType>& transport_mixture() const;

    // forward compatibility
    const TransportMixture<ThermoEvaluator,CoeffType>& mixture() const;

    // backward compatibility
    const std::vector<Viscosity*> & species_viscosities() const;

    //! Formatted print, by default to \p std::cout
    void print(std::ostream& os = std::cout) const;

    //! Formatted print.
    friend std::ostream& operator<<(std::ostream& os, const MixtureViscosity& mu)
    {
      mu.print(os);
      return os;
    }

  protected:

    const TransportMixture<ThermoEvaluator,CoeffType>& _transport_mixture;

    std::vector<Viscosity*> _species_viscosities;

  };

  template<typename Viscosity, typename ThermoEvaluator, class CoeffType>
  MixtureViscosity<Viscosity,ThermoEvaluator,CoeffType>::MixtureViscosity(  const TransportMixture<ThermoEvaluator,CoeffType>& transport_mixture )
    :  _transport_mixture(transport_mixture),
       _species_viscosities( transport_mixture.n_species(), NULL )
  {
    antioch_deprecated();
    return;
  }

  template<typename Viscosity, typename ThermoEvaluator, class CoeffType>
  MixtureViscosity<Viscosity,ThermoEvaluator,CoeffType>::~MixtureViscosity()
  {
    // Need to delete all the species viscosities we allocated
    for( typename std::vector<Viscosity*>::iterator it = _species_viscosities.begin();
	 it != _species_viscosities.end(); ++it )
      {
	delete (*it);
      }
    return;
  }

  template<typename Viscosity, typename ThermoEvaluator, class CoeffType>
  void MixtureViscosity<Viscosity,ThermoEvaluator,CoeffType>::add( const std::string& species_name,
						   const std::vector<CoeffType>& coeffs )
  {
    antioch_assert( _transport_mixture.species_name_map().find(species_name) !=
		    _transport_mixture.species_name_map().end() );

    unsigned int s = _transport_mixture.species_name_map().find(species_name)->second;

    antioch_assert_less_equal( s, _species_viscosities.size() );
    antioch_assert( !_species_viscosities[s] );

    _species_viscosities[s] = new Viscosity(coeffs);
    return;
  }

  template<typename Viscosity, typename ThermoEvaluator, class CoeffType>
  void MixtureViscosity<Viscosity,ThermoEvaluator,CoeffType>::reset_coeffs( const unsigned int s,
                                                            const std::vector<CoeffType> coeffs )
  {
    _species_viscosities[s]->reset_coeffs(coeffs);
  }

  template<typename Viscosity, typename ThermoEvaluator, class CoeffType>
  template<typename StateType>
  inline
  StateType MixtureViscosity<Viscosity,ThermoEvaluator,CoeffType>::operator()( const unsigned int s,
							       const StateType& T ) const
  {
    antioch_assert_less_equal( s, _species_viscosities.size() );
    antioch_assert( _species_viscosities[s] );

    return (*_species_viscosities[s])(T);
  }

  template<typename Viscosity, typename ThermoEvaluator, class CoeffType>
  inline
  const ChemicalMixture<CoeffType>& MixtureViscosity<Viscosity,ThermoEvaluator,CoeffType>::chemical_mixture() const
  {
    return _transport_mixture.chemical_mixture();
  }

  template<typename Viscosity, typename ThermoEvaluator, class CoeffType>
  inline
  const TransportMixture<ThermoEvaluator,CoeffType>& MixtureViscosity<Viscosity,ThermoEvaluator,CoeffType>::transport_mixture() const
  {
    return _transport_mixture;
  }

  template<typename Viscosity, typename ThermoEvaluator, class CoeffType>
  inline
  const TransportMixture<ThermoEvaluator,CoeffType>& MixtureViscosity<Viscosity,ThermoEvaluator,CoeffType>::mixture() const
  {
    return _transport_mixture;
  }

  template<typename Viscosity, typename ThermoEvaluator, class CoeffType>
  inline
  const std::vector<Viscosity*>& MixtureViscosity<Viscosity,ThermoEvaluator,CoeffType>::species_viscosities() const
  {
    return _species_viscosities;
  }

  template<typename Viscosity, typename ThermoEvaluator, class CoeffType>
  void MixtureViscosity<Viscosity,ThermoEvaluator,CoeffType>::print(std::ostream& os) const
  {
    antioch_assert_equal_to( _species_viscosities.size(), _transport_mixture.n_species() );

    for( unsigned int s = 0; s < _transport_mixture.n_species(); s++ )
      {
	const Species& species = _transport_mixture.species_list()[s];
	const std::string& name = _transport_mixture.species_inverse_name_map().find( species  )->second;

	os << "mu(" << name << ") = " << (*_species_viscosities[s]) << std::endl;
      }

    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_MIXTURE_VISCOSITY_H
