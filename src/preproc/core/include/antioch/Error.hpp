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

#ifndef _ERROR_HANDLING_
#define _ERROR_HANDLING_

//
#include <string>
#include <iostream>
#include <cstdlib>

namespace Antioch
{
/*! \class Error
 * \brief class Error, to deal with errors and warning, to make it simple
 *
 * The class that takes the errors. Only a printing option at this
 * time. Takes the type of error (warning, error, whatever),
 * the name of the method and the message, all strings.
 */

class Error{
  public:
/*! \brief Default constructor, does nothing, useless*/
    Error(){exit(1);}
/*! \brief Constructor that should be used
 *
 * \param const string &type, error, warning, note of death,...
 * \param const string &method, name of the method returning the error
 * \param const string &message, the message to be printed
 * \return */
    Error(const std::string &type,const std::string &method, const std::string &message);
/*! \brief Default destructor*/
    ~Error(){}

  private:
/*! \brief Prints the error on the screen, using the constructor's parameters*/
    void printError(const std::string &type,const std::string &method,const std::string &message) const;
/*! \brief Sorting between error and warning*/
    bool isError(const std::string &type) const;

};

/*!\brief Small public method to shorten warning message generation.
 *
 * It creates an Error instance for warning.
 */
void antiochWarning(const std::string &method,const std::string &message);
/*!\brief Small public method to shorten error message generation and stops program.
 *
 * It creates an Error instance for error.
 */
void antiochError(const std::string &method,const std::string &message);
}

#endif
