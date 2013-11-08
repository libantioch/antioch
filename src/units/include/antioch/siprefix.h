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
#ifndef ANTIOCH_SI_PREFIX_H
#define ANTIOCH_SI_PREFIX_H

//Antioch
#include "antioch/metaprogramming.h"
#include "antioch/antioch_asserts.h"

//C++
#include <string>
#include <iostream>

namespace Antioch{

/*!\class SIPrefixes
 * \brief Prefixes in unit
 *
 * This class associates a std::string and
 * a double. We store here the prefixes and the
 * associated value defined in the file unit_defs.hpp
 *
 * To add prefixes, one needs to add a SIPrefixes instance in the
 * const SIPrefixes Prefixes[] variable in the unit_defs.hpp file:
 *
 * SIPrefixes("prefix",value)
 */
template <typename T = double>
class SIPrefixes{
    public:
/*! \brief Copy constructor, uses the SIPrefixes &operator=(const SIPrefixe&)*/
      template <typename P>
      SIPrefixes(const SIPrefixes<P> &rhs){*this = rhs;}
/*! \brief Default constructor, a std::string and a value*/
      SIPrefixes(const std::string &str, const std::string & na, const T &num):
                _symbol(str),_name(na),_value(num){}
/*! \brief Default destructor*/
      ~SIPrefixes(){}

/*! \brief Value getter*/
      template <typename P = T>
      const P            value() const {return _value;}
/*! \brief Symbol getter*/
      const std::string symbol() const {return _symbol;}
/*! \brief Name getter*/
      const std::string   name() const {return _name;}

/*! \brief Assignement operator, copy value and std::string*/
      template <typename P>
      SIPrefixes &operator=(const SIPrefixes<P> &rhs);
/*! \brief Comparison operator, equal if values are equal*/
      template <typename P>
      bool const operator==(const SIPrefixes<P> &rhs) const {return (_value == rhs.value());}
/*! \brief Comparison operator, not equal is not "equal"*/
      template <typename P>
      bool const operator!=(const SIPrefixes<P> &rhs) const {return !(*this == rhs);}

    private:
/*! \brief Default constructor, never uninitialize this stuff*/
      SIPrefixes(){antioch_error();}
/*! \brief Two std::strings for the symbol and the name*/
      std::string _symbol,_name;
/*! \brief A double for the value*/
      T _value;
};

template <typename T>
template <typename P>
inline
SIPrefixes<T> &SIPrefixes<T>::operator=(const SIPrefixes<P> &rhs)
{
  if(this == &rhs){return *this;}
  _symbol = rhs.symbol();
  Antioch::init_clone(_value,rhs.value());
  return *this;
}
}// Antioch namespace

#endif
