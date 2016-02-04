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

#include "antioch/string_utils.h"

namespace Antioch
{
  void split_string( const std::string& input,
                     const std::string& delimiter,
                     std::vector<std::string>& results )
  {
    // Skip delimiters at beginning.
    std::string::size_type first_pos = input.find_first_not_of(delimiter, 0);

    std::string::size_type pos = input.find_first_of(delimiter, first_pos);

    while (std::string::npos != pos || std::string::npos != first_pos)
      {
        // Found a token, add it to the vector.
        results.push_back(input.substr(first_pos, pos - first_pos));

        // Skip delimiters.  Note the "not_of"
        first_pos = input.find_first_not_of(delimiter, pos);

        // Find next delimiter
        pos = input.find_first_of(delimiter, first_pos);
      }
  }

}
