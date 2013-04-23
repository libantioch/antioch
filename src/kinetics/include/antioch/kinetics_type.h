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

#ifndef _ANTIOCH_KINETICS_TYPE_H
#define _ANTIOCH_KINETICS_TYPE_H

//Antioch
#include "antioch/kinetics_enum.h"

//C++
#include <string>
#include <iostream>

namespace Antioch{
/*!
 *
 * \class KineticsType
 * \brief base class for kinetics models
 */
template <typename CoeffType>
class KineticsType{
   public:
      KineticsType(const KinMod::KinMod type);
      virtual ~KineticsType();

    //!\return error because I cannot make it pure virtual.
    template <typename StateType>
    StateType operator()(const StateType& T) const;

    //!\return error because I cannot make it pure virtual.
    template <typename StateType>
    StateType derivative( const StateType& T ) const;

    //!\return error because I cannot make it pure virtual.
    template <typename StateType>
    void rate_and_derivative(const StateType& T, StateType& rate, StateType& drate_dT) const;

    virtual const std::string numeric() const = 0;

    //! Formatted print, by default to \p std::cout
    void print(std::ostream& os = std::cout) const;

    //! Formatted print.
    friend std::ostream& operator<<(std::ostream& os, const KineticsType& rate)
    {
      rate.print(os);
      return os;
    }

   private:
   KinMod::KinMod my_type;
};

  /* ------------------------- Inline Functions -------------------------*/
  template <typename CoeffType>
  inline
  void KineticsType<CoeffType>::print(std::ostream& os) const
  {
    os << numeric() << std::endl;
  }

  template <typename CoeffType>
  inline
  KineticsType<CoeffType>::KineticsType(const KinMod::KinMod type):
    my_type(type)
  {
      return;
  }

  template <typename CoeffType>
  inline
  KineticsType<CoeffType>::~KineticsType()
  {
     return;
  }
}

#endif
