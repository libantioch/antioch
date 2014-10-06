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

#ifndef ANTIOCH_TRANSPORT_MIXTURE_H
#define ANTIOCH_TRANSPORT_MIXTURE_H

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/transport_species.h"
#include "antioch/metaprogramming.h"
#include "antioch/chemical_mixture.h"
#include "antioch/transport_species_ascii_parsing.h"
#include "antioch/default_filename.h"

// C++
#include <vector>
#include <string>

namespace Antioch
{
  //! Class storing chemical mixture properties
  /*!
    This class manages the list of TransportSpecies for a requested set
    of species from input.
  */
  template<typename ThermoEvaluator,typename CoeffType=double>
  class TransportMixture
  {

     typedef unsigned int Species;

  public:
    
    TransportMixture( const ChemicalMixture<CoeffType> &mixture, const ThermoEvaluator & t,
                      const std::string & filename = DefaultFilename::transport_mixture());

    ~TransportMixture();

    //! ChemicalMixture method
    unsigned int n_species() const;

    //! ChemicalMixture method
    const std::vector<Species>& species_list() const;

    //! ChemicalMixture method
    const std::map<Species,std::string>& species_inverse_name_map() const;

    //! ChemicalMixture method
    const std::map<std::string, Species>& species_name_map() const;

    //! \returns the chemical mixture
    const ChemicalMixture<CoeffType> & chemical_mixture() const;

    //! \returns the thermodynamics evaluator
    const ThermoEvaluator & thermo() const;

    void add_species( const unsigned int index,
                      CoeffType LJ_depth, CoeffType LJ_diameter,
                      CoeffType dipole_moment, CoeffType polarizability, CoeffType rotational_relaxation, CoeffType mass);

    const std::vector<TransportSpecies<CoeffType>*>& transport_species() const;

  protected:

    const ChemicalMixture<CoeffType>        & _chemical_mixture;

    const ThermoEvaluator                   & _thermo;

    std::vector<TransportSpecies<CoeffType>*> _transport_species;

  private:
    TransportMixture();

  };


  /* ------------------------- Inline Functions -------------------------*/

  template<typename ThermoEvaluator,typename CoeffType>
  inline
  unsigned int TransportMixture<ThermoEvaluator,CoeffType>::n_species() const
  {
    return _chemical_mixture.n_species();
  }

  template<typename ThermoEvaluator,typename CoeffType>
  inline
  const std::vector<Species>& TransportMixture<ThermoEvaluator,CoeffType>::species_list() const
  { 
    return _chemical_mixture.species_list();
  }

  template<typename ThermoEvaluator,typename CoeffType>
  inline
  const std::map<Species,std::string>& TransportMixture<ThermoEvaluator,CoeffType>::species_inverse_name_map() const
  {
    return _chemical_mixture.species_inverse_name_map();
  }

  template<typename ThermoEvaluator,typename CoeffType>
  inline
  const std::map<std::string,Species>& TransportMixture<ThermoEvaluator,CoeffType>::species_name_map() const
  {
    return _chemical_mixture.species_name_map();
  }

  template<typename ThermoEvaluator,typename CoeffType>
  inline
  const ChemicalMixture<CoeffType> & TransportMixture<ThermoEvaluator,CoeffType>::chemical_mixture() const
  {
     return _chemical_mixture;
  }

  template<typename ThermoEvaluator,typename CoeffType>
  inline
  const ThermoEvaluator & TransportMixture<ThermoEvaluator,CoeffType>::thermo() const
  {
     return _thermo;
  }

  template<typename ThermoEvaluator,typename CoeffType>
  inline
  const std::vector<TransportSpecies<CoeffType>*>& TransportMixture<ThermoEvaluator,CoeffType>::transport_species() const
  {
    return _transport_species;
  }

  template<typename ThermoEvaluator,typename CoeffType>
  inline
  TransportMixture<ThermoEvaluator,CoeffType>::TransportMixture( const ChemicalMixture<CoeffType>& chem_mix, const ThermoEvaluator & t, const std::string & filename )
    : _chemical_mixture( chem_mix),
      _thermo(t),
      _transport_species(_chemical_mixture.n_species(), NULL )
  {

    // Now read in transport properties for the requested species and stash
    read_transport_species_data_ascii(*this, filename);

    // check we have everyone requested
    for( unsigned int s = 0; s < _transport_species.size(); ++s )
      {
        if(!_transport_species[s]) // it is not mandatory, Blottner or Sutherland are self-sufficient
        {
           std::cerr << "Warning: missing transport data for species " << _chemical_mixture.species_inverse_name_map().at(
                                                                          _chemical_mixture.species_list()[s]) << "\n"
                     << "Be sure to use a transport model that does not require the default data"
                     << std::endl;
        }
      }

    return;

  }

  template<typename ThermoEvaluator,typename CoeffType>
  inline
  TransportMixture<ThermoEvaluator,CoeffType>::~TransportMixture()
  {
    // Clean up all the TransportSpecies we stored
    for( typename std::vector<TransportSpecies<CoeffType>* >::iterator it = _transport_species.begin();
         it < _transport_species.end(); ++it )
      {
        delete (*it);
      }

    return;
  }

  template<typename ThermoEvaluator,typename CoeffType>
  inline
  void TransportMixture<ThermoEvaluator,CoeffType>::add_species( const unsigned int index,
                                                CoeffType LJ_depth,
                                                CoeffType LJ_diameter,
                                                CoeffType dipole_moment,
                                                CoeffType polarizability, 
                                                CoeffType rotational_relaxation,
                                                CoeffType mass)
  {
    Species name_enum = _chemical_mixture.species_list()[index];
    _transport_species[index] =
      new TransportSpecies<CoeffType>(name_enum, LJ_depth, LJ_diameter, dipole_moment, polarizability, rotational_relaxation,mass);

    return;
  }

} // end namespace Antioch

#endif //ANTIOCH_TRANSPORT_MIXTURE_H
