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
#ifndef _SI_PREFIX_
#define _SI_PREFIX_

//Antioch

//C++
#include <string>

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
 * SIPrefixes("prefixe",value)
 */
class SIPrefixes{
    public:
/*! \brief Copy constructor, uses the SIPrefixes &operator=(const SIPrefixe&)*/
      SIPrefixes(const SIPrefixes &rhs){*this = rhs;}
/*! \brief Default constructor, a std::string and a value*/
      SIPrefixes(std::string str="",std::string na="",double num=0.):_symbol(str),_name(na),_value(num){}
/*! \brief Default destructor*/
      ~SIPrefixes(){}

/*! \brief Value getter*/
      double value()             const {return _value;}
/*! \brief Symbol getter*/
      const std::string symbol() const {return _symbol;}
/*! \brief Name getter*/
      const std::string name()   const {return _name;}

/*! \brief Assignement operator, copy value and std::string*/
      SIPrefixes &operator=(const SIPrefixes &rhs);
/*! \brief Comparison operator, equal if values are equal*/
      bool const operator==(const SIPrefixes &rhs) const {return (_value == rhs.value());}
/*! \brief Comparison operator, not equal is not "equal"*/
      bool const operator!=(const SIPrefixes &rhs) const {return !(*this == rhs);}

    private:
/*! \brief Two std::strings for the symbol and the name*/
      std::string _symbol,_name;
/*! \brief A double for the value*/
      double _value;
};

}

#endif
