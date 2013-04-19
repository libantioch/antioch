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

namespace Antioch{
/*!
 *
 * \class KineticsType
 * \brief base class for kinetics models
 */
class KineticsType{
   public:
      KineticsType();
      virtual ~KineticsType();

    //!\return error because I cannot make it pure virtual.
    template <typename StateType>
    virtual StateType operator()(const StateType& T) const;

    //!\return error because I cannot make it pure virtual.
    template <typename StateType>
    virtual StateType derivative( const StateType& T ) const;

    //!\return error because I cannot make it pure virtual.
    template <typename StateType>
    virtual void rate_and_derivative(const StateType& T, StateType& rate, StateType& drate_dT) const;
};

  /* ------------------------- Inline Functions -------------------------*/
  inline
  KineticsType::KineticsType()
  {
      return;
  }

  inline
  KineticsType::~KineticsType()
  {
     return;
  }

  inline
  template <typename StateType>
  StateType operator()(const StateType& T) const
  {
    std::cerr << "Kinetics base class, get outta here!!" << std::endl;
    antioch_error();
  }

  inline
  template <typename StateType>
  StateType derivative( const StateType& T ) const
  {
    std::cerr << "Kinetics base class, get outta here!!" << std::endl;
    antioch_error();
  }

  inline
  template <typename StateType>
  void rate_and_derivative(const StateType& T, StateType& rate, StateType& drate_dT) const
  {
    std::cerr << "Kinetics base class, get outta here!!" << std::endl;
    antioch_error();
  }

}

#endif
