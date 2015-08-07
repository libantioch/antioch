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

#include "antioch/transport_species_parsing.h"

// Antioch
#include "antioch/parser_base.h"

namespace Antioch
{
  template <typename NumericType>
  void read_transport_species_data(ParserBase<NumericType> * parser, TransportMixture<NumericType> & transport)
  {
    parser->read_transport_data(transport);

    // sanity check, we may require these informations
    bool fail(false);
    for(unsigned int s = 0; s < transport.n_species(); s++)
      {
        if(!transport.transport_species()[s])
          {
            fail = true;
            break;
          }
      }
    if(fail)
      {
        std::cerr << "Molecule(s) is(are) missing in transport description.  Please update the information."
                  << "  Currently using file " << parser->file() << ".\n"
                  << "You might have some problem later if you need these description.  "
                  << "Missing molecule(s) is(are):" << std::endl;
        for(unsigned int i = 0; i < transport.n_species(); i++)
          {
            if(!transport.transport_species()[i])
              {
                std::cerr << transport.species_inverse_name_map().at(i) << std::endl;
              }
          }
      }
  }

  // Instantiate
  template void read_transport_species_data<float>( ParserBase<float>*, TransportMixture<float>& );
  template void read_transport_species_data<double>( ParserBase<double>*, TransportMixture<double>& );
  template void read_transport_species_data<long double>( ParserBase<long double>*, TransportMixture<long double>& );

} // end namespace Antioch
