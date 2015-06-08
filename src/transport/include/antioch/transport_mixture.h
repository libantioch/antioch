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
#include "antioch/chemical_mixture.h"
#include "antioch/parsing_enum.h"

// C++
#include <vector>
#include <string>
#include <map>

namespace Antioch
{
  //! Class storing chemical mixture properties
  /*!
    This class manages the list of TransportSpecies for a requested set
    of species from input.
  */
  template<typename CoeffType=double>
  class TransportMixture
  {

     typedef unsigned int Species;

  public:

    TransportMixture( const ChemicalMixture<CoeffType> &mixture,
                      const std::string & filename = DefaultFilename::transport_mixture(),
                      bool verbose = true, ParsingType type = ASCII );

    TransportMixture( const ChemicalMixture<CoeffType> &mixture,
                      ParserBase<CoeffType> * parser);

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

    void add_species( const unsigned int index,
                      CoeffType LJ_depth, CoeffType LJ_diameter,
                      CoeffType dipole_moment, CoeffType polarizability, CoeffType rotational_relaxation, CoeffType mass);

    const std::vector<TransportSpecies<CoeffType>*>& transport_species() const;

    const TransportSpecies<CoeffType>& transport_species( unsigned int s ) const;

  protected:

    const ChemicalMixture<CoeffType>& _chemical_mixture;

    std::vector<TransportSpecies<CoeffType>*> _transport_species;

  private:
    TransportMixture();

  };


  /* ------------------------- Inline Functions -------------------------*/

  template<typename CoeffType>
  inline
  unsigned int TransportMixture<CoeffType>::n_species() const
  {
    return _chemical_mixture.n_species();
  }

  template<typename CoeffType>
  inline
  const std::vector<Species>& TransportMixture<CoeffType>::species_list() const
  {
    return _chemical_mixture.species_list();
  }

  template<typename CoeffType>
  inline
  const std::map<Species,std::string>& TransportMixture<CoeffType>::species_inverse_name_map() const
  {
    return _chemical_mixture.species_inverse_name_map();
  }

  template<typename CoeffType>
  inline
  const std::map<std::string,Species>& TransportMixture<CoeffType>::species_name_map() const
  {
    return _chemical_mixture.species_name_map();
  }

  template<typename CoeffType>
  inline
  const ChemicalMixture<CoeffType> & TransportMixture<CoeffType>::chemical_mixture() const
  {
     return _chemical_mixture;
  }

  template<typename CoeffType>
  inline
  const std::vector<TransportSpecies<CoeffType>*>& TransportMixture<CoeffType>::transport_species() const
  {
    return _transport_species;
  }

  template<typename CoeffType>
  inline
  const TransportSpecies<CoeffType>& TransportMixture<CoeffType>::transport_species(unsigned int s) const
  {
    antioch_assert_less(s,_transport_species.size());
    antioch_assert(_transport_species[s]);

    return (*_transport_species[s]);
  }

} // end namespace Antioch

#endif //ANTIOCH_TRANSPORT_MIXTURE_H
