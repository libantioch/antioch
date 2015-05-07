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
#include "antioch/parsing_enum.h"
#include "antioch/transport_species_parsing.h"
#include "antioch/default_filename.h"
#include "antioch/parser_base.h"

// C++
#include <vector>
#include <string>

namespace Antioch
{

  template <typename CoeffType>
  class ASCIIParser;

  template <typename CoeffType>
  class ChemKinParser;

  template <typename CoeffType>
  class XMLParser;

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

  protected:

    const ChemicalMixture<CoeffType>        & _chemical_mixture;

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
  TransportMixture<CoeffType>::TransportMixture( const ChemicalMixture<CoeffType>& chem_mix,
                                                                 const std::string & filename, bool verbose, ParsingType type )
    : _chemical_mixture( chem_mix),
      _transport_species(_chemical_mixture.n_species(), NULL )
  {

   ParserBase<CoeffType> * parser(NULL);
    switch(type)
    {
      case ASCII:
         parser = new ASCIIParser<CoeffType>(filename,verbose);
         break;
      case CHEMKIN:
         parser = new ChemKinParser<CoeffType>(filename,verbose);
         break;
      case XML:
         parser = new XMLParser<CoeffType>(filename,verbose);
         break;
      default:
         antioch_parsing_error("unknown type");
    }

    // Now read in transport properties for the requested species and stash
    read_transport_species_data(parser,*this);

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

    delete parser;

    return;

  }

  template<typename CoeffType>
  inline
  TransportMixture<CoeffType>::TransportMixture( const ChemicalMixture<CoeffType> &mixture,
                                                                 ParserBase<CoeffType> * parser)
    : _chemical_mixture( mixture),
      _transport_species(_chemical_mixture.n_species(), NULL )
  {
    // Now read in transport properties for the requested species and stash
    read_transport_species_data(parser,*this);

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

  template<typename CoeffType>
  inline
  TransportMixture<CoeffType>::~TransportMixture()
  {
    // Clean up all the TransportSpecies we stored
    for( typename std::vector<TransportSpecies<CoeffType>* >::iterator it = _transport_species.begin();
         it < _transport_species.end(); ++it )
      {
        delete (*it);
      }

    return;
  }

  template<typename CoeffType>
  inline
  void TransportMixture<CoeffType>::add_species( const unsigned int index,
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
