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
//
#ifndef ANTIOCH_UNIT_BASE_H
#define ANTIOCH_UNIT_BASE_H

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/insi.h"
#include "antioch/converter.h"

//C++


namespace Antioch{

//!wrapper for unit storage
template <typename T = double>
class UnitBase{
   public:
     UnitBase(const std::string &sym, const std::string &nam, // name & symbol
              const T &fac,   const T &trans, // converter
              int mi,int kgi=0, int si=0, int Ai=0, int Ki=0, int moli=0, int cdi=0, int radi=0):// power array
                _conversion(fac,trans),
                _power_base(mi,kgi,si,Ai,Ki,moli,cdi,radi),
                _symbol(sym),_name(nam){}

     ~UnitBase(){}


    template <typename P = T>
    const P               factor() const {return _conversion.geta();}
    template <typename P = T>
    const P           translator() const {return _conversion.getb();}
    template <typename P = T>
    const Converter<P> converter() const {return _conversion;}
    const std::string       name() const {return _name;}
    const std::string     symbol() const {return _symbol;}
    const InSI       power_array() const {return _power_base;}

   private:
//!no default possible
     UnitBase(){antioch_error();}

     Converter<T> _conversion;
     InSI         _power_base;
     std::string  _symbol;
     std::string  _name;
};

}// end namespace
#endif
