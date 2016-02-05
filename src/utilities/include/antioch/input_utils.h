//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014-2016 Paul T. Bauman, Benjamin S. Kirk,
//                         Sylvain Plessis, Roy H. Stonger
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

#ifndef ANTIOCH_INPUT_UTILS_H
#define ANTIOCH_INPUT_UTILS_H

#include <istream>
#include <locale> // isblank

namespace Antioch
{
  /*!
    Skip comment lines in the header of an ASCII
    text file prefixed with the comment character
    'comment_start'. This is put in Antioch namespace so we
    can reuse in a few classes where we are reading in text
    tables.
    Originally taken from FIN-S.
  */
  void skip_comment_lines( std::istream &in, const char comment_start);

  inline
  void skip_comment_lines( std::istream &in, const char comment_start)
  {
    char c, line[256];

    in.get(c);
    while(std::isblank(c))in.get(c);

    in.putback(c);

    while (in.get(c), c==comment_start)
       in.getline (line, 255);

    // put back first character of
    // first non-comment line
    in.putback (c);
  }

} // end namespace Antioch

#endif //ANTIOCH_INPUT_UTILS_H
