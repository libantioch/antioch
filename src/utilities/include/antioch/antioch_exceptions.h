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

#ifndef ANTIOCH_EXCEPTIONS_H
#define ANTIOCH_EXCEPTIONS_H

// C++
#include <stdexcept>
#include <string>

namespace Antioch
{
  /*!
   * A class to represent the internal "this should never happen"
   * errors, to be thrown by "antioch_error();"
   */
  class LogicError : public std::logic_error
  {
  public:
    LogicError() : std::logic_error( "Error in Antioch internal logic" ) {}
  };

  /*!
   * A class to stub for features that should be in Antioch, but
   * haven't been written yet, to be thrown by
   * "antioch_not_implemented();"
   */
  class NotImplemented : public std::logic_error
  {
    public:
    NotImplemented() : std::logic_error( "Error: not implemented!" ) {}
  };
  
  /*!
   * A class representing a failed attempt by the library to open a
   * file (or construct an fstream, etc), to be thrown by
   * "antioch_file_error(filename);" For ease of debugging, "filename"
   * should include any (absolute or relative or implicit) pathname
   * that was part of the failed open.
   */
  class FileError : public std::runtime_error
  {
  public:
    FileError(const std::string& filename) : std::runtime_error( "Error accessing file: " + filename ) {}
  };

} // end namespace Antioch

#define ANTIOCH_THROW(e) do { throw e; } while (0)

#endif // ANTIOCH_EXCEPTIONS_H
