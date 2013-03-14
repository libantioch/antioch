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

#ifndef ANTIOCH_MIXTURE_THERMAL_CONDUCTIVITY_H
#define ANTIOCH_MIXTURE_THERMAL_CONDUCTIVITY_H

// C++
#include <string>
#include <vector>

// Antioch
#include "antioch/chemical_mixture.h"

namespace Antioch
{
  template<typename ThermalConductivity, class CoeffType=double>
  class MixtureThermalConductivity
  {
  public:

    MixtureThermalConductivity( const ChemicalMixture<CoeffType>& chem_mixture );
    ~MixtureThermalConductivity();

    template <typename StateType>
    StateType operator()( const unsigned int s, const StateType T ) const;

    void add( const std::string& species_name,
	      const std::vector<CoeffType>& coeffs );

    const ChemicalMixture<CoeffType>& chemical_mixture() const;

    //! Formatted print, by default to \p std::cout
    void print(std::ostream& os = std::cout) const;

    //! Formatted print.
    friend std::ostream& operator<<(std::ostream& os, const MixtureThermalConductivity& k)
    {
      k.print(os);
      return os;
    }

  protected:

    const ChemicalMixture<CoeffType>& _chem_mixture;

    std::vector<ThermalConductivity*> _species_conductivities;

  };

  template<typename ThermalConductivity, class CoeffType>
  MixtureThermalConductivity<ThermalConductivity,CoeffType>::MixtureThermalConductivity( const ChemicalMixture<CoeffType>& chem_mixture )
    :  _chem_mixture(chem_mixture),
       _species_conductivities( chem_mixture.n_species(), NULL )
  {
    return;
  }

  template<typename ThermalConductivity, class CoeffType>
  MixtureThermalConductivity<ThermalConductivity,CoeffType>::~MixtureThermalConductivity()
  {
    // Need to delete all the species conductivities we allocated
    for( typename std::vector<ThermalConductivity*>::iterator it = _species_conductivities.begin();
	 it != _species_conductivities.end(); ++it )
      {
	delete (*it);
      }
    return;
  }

  template<typename ThermalConductivity, class CoeffType>
  void MixtureThermalConductivity<ThermalConductivity,CoeffType>::add( const std::string& species_name,
						   const std::vector<CoeffType>& coeffs )
  {
    antioch_assert( _chem_mixture.active_species_name_map().find(species_name) !=
		    _chem_mixture.active_species_name_map().end() );

    unsigned int s = _chem_mixture.active_species_name_map().find(species_name)->second;

    antioch_assert_less_equal( s, _species_conductivities.size() );
    antioch_assert( !_species_conductivities[s] );

    _species_conductivities[s] = new ThermalConductivity(coeffs);

    return;
  }

  template<typename ThermalConductivity, class CoeffType>
  template<typename StateType>
  inline
  StateType MixtureThermalConductivity<ThermalConductivity,CoeffType>::operator()( const unsigned int s,
										   const StateType T ) const
  {
    antioch_assert_less_equal( s, _species_conductivities.size() );
    antioch_assert( _species_conductivities[s] );

    return (*_species_conductivities[s])(T);
  }

  template<typename ThermalConductivity, class CoeffType>
  inline
  const ChemicalMixture<CoeffType>& MixtureThermalConductivity<ThermalConductivity,CoeffType>::chemical_mixture() const
  {
    return _chem_mixture;
  }

  template<typename ThermalConductivity, class CoeffType>
  void MixtureThermalConductivity<ThermalConductivity,CoeffType>::print(std::ostream& os) const
  {
    antioch_assert_equal_to( _species_conductivities.size(), _chem_mixture.n_species() );

    for( unsigned int s = 0; s < _chem_mixture.n_species(); s++ )
      {
	const Species& species = _chem_mixture.species_list()[s];
	const std::string& name = _chem_mixture.species_inverse_name_map().find( species  )->second;

	os << "k(" << name << ") = " << (*_species_conductivities[s]) << std::endl;
      }

    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_MIXTURE_THERMAL_CONDUCTIVITY_H
