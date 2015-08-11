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

// This class
#include "antioch/parser_base.h"

// Antioch
#include "antioch/antioch_numeric_type_instantiate_macro.h"

namespace Antioch
{
  template <typename NumericType>
  ParserBase<NumericType>::ParserBase(const std::string & type, const std::string & file, bool verbose, const std::string & comments)
    : _type(type),
      _file(file),
      _verbose(verbose),
      _comments(comments)
  {
    std::stringstream os;
    os <<  "\n*********************************************************\n"
       << "This method is not available with a " << _type << " parser.\n"
       << "Parsing file " << _file << ".\n"
       << "No format has been defined yet.  Maybe contribute?\n"
       << "https://github.com/libantioch/antioch\n"
       << "\n\n*********************************************************\n\n";

    _not_implemented = os.str();
  }

  template <typename NumericType>
  void ParserBase<NumericType>::skip_comments(std::istream & doc)
  {
    for(unsigned int c = 0; c < _comments.size(); c++)
      {
        skip_comment_lines(doc, _comments[c]);
      }
  }

  template <typename NumericType>
  ParsingType ParserBase<NumericType>::enum_type() const
  {
    ParsingType PType(ASCII);
    if(_type == "ascii")
      {
        PType = ASCII;
      }
    else if(_type == "ChemKin")
      {
        PType = CHEMKIN;
      }
    else if(_type == "XML")
      {
        PType = XML;
      }
    else
      {
        antioch_parsing_error(std::string("unknown parser type!!! " + _type));
      }

    return PType;
  }

} // end namespace Antioch

// Instantiate
ANTIOCH_NUMERIC_TYPE_CLASS_INSTANTIATE(Antioch::ParserBase);
