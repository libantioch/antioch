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

#ifndef ANTIOCH_TRANSPORT_ASCII_PARSING_H
#define ANTIOCH_TRANSPORT_ASCII_PARSING_H

// Antioch
#include "antioch/input_utils.h"
#include "antioch/physical_constants.h"

//C++
#include <iostream>
#include <string>
#include <fstream>


namespace Antioch{

  //Forward declaration
  template <typename ThermoEvaluator, typename NumericType>
  class TransportMixture;


  template <typename ThermoEvaluator, typename NumericType>
  void read_transport_species_data_ascii(TransportMixture<ThermoEvaluator,NumericType> & transport, const std::string & filename);

/*----------- inline functions ----------------*/


  template <typename ThermoEvaluator, typename NumericType>
  void read_transport_species_data_ascii(TransportMixture<ThermoEvaluator,NumericType> & transport, const std::string & filename)
  {
    std::ifstream in(filename.c_str());
    if(!in.is_open())
    {
      std::cerr << "ERROR: unable to load file " << filename << std::endl;
      antioch_error();
    }


    skip_comment_lines(in, '#');

    std::string name;
    NumericType LJ_eps_kB;
    NumericType LJ_sigma;
    NumericType dipole_moment;
    NumericType pol;
    NumericType Zrot;

    while (in.good())
      {
          in >> name >> LJ_eps_kB >> LJ_sigma >> dipole_moment>> pol >> Zrot;
          if(transport.chemical_mixture().species_name_map().count(name))
          {
              unsigned int place = transport.chemical_mixture().species_name_map().at(name);
                // TODO: better unit checking
              NumericType mass = transport.chemical_mixture().M(place) * NumericType(1e-3); //to SI
// adding species in mixture
              transport.add_species(place,LJ_eps_kB,LJ_sigma,dipole_moment,pol,Zrot,mass);
          }
      }
     in.close();
    return;
  }

} //end namespace Antioch

#endif
