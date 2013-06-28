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
#include "antioch/antioch_asserts.h"

//C++
#include <string>
#include <iostream>

namespace Antioch{

template <typename CoeffType>
class HercourtEssenRate;
template <typename CoeffType>
class BerthelotRate;
template <typename CoeffType>
class ArrheniusRate;
template <typename CoeffType>
class BerthelotHercourtEssenRate;
template <typename CoeffType>
class KooijRate;
template <typename CoeffType>
class VantHoffRate;

/*!
 *
 * \class KineticsType
 * \brief base class for kinetics models
 */
template <typename CoeffType>
class KineticsType{
   public:
      KineticsType(const KineticsModel::KineticsModel type);
      virtual ~KineticsType();

    //!\return error because I cannot make it pure virtual.
    template <typename StateType>
    StateType operator()(const StateType& T) const;

    //!\return error because I cannot make it pure virtual.
    template <typename StateType>
    StateType derivative( const StateType& T ) const;

    //!\return error because I cannot make it pure virtual.
    template <typename StateType>
    void compute_rate_and_derivative(const StateType& T, StateType& rate, StateType& drate_dT) const;

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
   KineticsModel::KineticsModel my_type;
};

  /* ------------------------- Inline Functions -------------------------*/
  template <typename CoeffType>
  inline
  void KineticsType<CoeffType>::print(std::ostream& os) const
  {
    os << numeric();
  }

  template <typename CoeffType>
  inline
  KineticsType<CoeffType>::KineticsType(const KineticsModel::KineticsModel type):
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

  template <typename CoeffType>
  template <typename StateType>
  inline
  StateType KineticsType<CoeffType>::operator()(const StateType& T) const
  {
      switch (my_type) 
      {
      case KineticsModel::HERCOURT_ESSEN:
        return (static_cast<const HercourtEssenRate<CoeffType>*>(this))->HercourtEssenRate<CoeffType>::rate(T);
        break;
      case KineticsModel::BERTHELOT:
        return (static_cast<const BerthelotRate<CoeffType>*>(this))->BerthelotRate<CoeffType>::rate(T);
        break;
      case KineticsModel::ARRHENIUS:
        return (static_cast<const ArrheniusRate<CoeffType>*>(this))->ArrheniusRate<CoeffType>::rate(T);
        break;
      case KineticsModel::BHE:
        return (static_cast<const BerthelotHercourtEssenRate<CoeffType>*>(this))->BerthelotHercourtEssenRate<CoeffType>::rate(T);
        break;
      case KineticsModel::KOOIJ:
        return (static_cast<const KooijRate<CoeffType>*>(this))->KooijRate<CoeffType>::rate(T);
        break;
      case KineticsModel::VANTHOFF:
        return (static_cast<const VantHoffRate<CoeffType>*>(this))->VantHoffRate<CoeffType>::rate(T);
        break;
      default:
        antioch_error();
        break;
      }
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  StateType KineticsType<CoeffType>::derivative( const StateType& T ) const
  {
      switch (my_type) 
      {
      case KineticsModel::HERCOURT_ESSEN:
        return (static_cast<const HercourtEssenRate<CoeffType>*>(this))->HercourtEssenRate<CoeffType>::derivative(T);
        break;
      case KineticsModel::BERTHELOT:
        return (static_cast<const BerthelotRate<CoeffType>*>(this))->BerthelotRate<CoeffType>::derivative(T);
        break;
      case KineticsModel::ARRHENIUS:
        return (static_cast<const ArrheniusRate<CoeffType>*>(this))->ArrheniusRate<CoeffType>::derivative(T);
        break;
      case KineticsModel::BHE:
        return (static_cast<const BerthelotHercourtEssenRate<CoeffType>*>(this))->BerthelotHercourtEssenRate<CoeffType>::derivative(T);
        break;
      case KineticsModel::KOOIJ:
        return (static_cast<const KooijRate<CoeffType>*>(this))->KooijRate<CoeffType>::derivative(T);
        break;
      case KineticsModel::VANTHOFF:
        return (static_cast<const VantHoffRate<CoeffType>*>(this))->VantHoffRate<CoeffType>::derivative(T);
        break;
      default:
        antioch_error();
        break;
      }
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void KineticsType<CoeffType>::compute_rate_and_derivative(const StateType& T, StateType& rate, StateType& drate_dT) const
  {
      switch (my_type) 
      {
      case KineticsModel::HERCOURT_ESSEN:
        return (static_cast<const HercourtEssenRate<CoeffType>*>(this))->HercourtEssenRate<CoeffType>::rate_and_derivative(T,rate,drate_dT);
        break;
      case KineticsModel::BERTHELOT:
        return (static_cast<const BerthelotRate<CoeffType>*>(this))->BerthelotRate<CoeffType>::rate_and_derivative(T,rate,drate_dT);
        break;
      case KineticsModel::ARRHENIUS:
        return (static_cast<const ArrheniusRate<CoeffType>*>(this))->ArrheniusRate<CoeffType>::rate_and_derivative(T,rate,drate_dT);
        break;
      case KineticsModel::BHE:
        return (static_cast<const BerthelotHercourtEssenRate<CoeffType>*>(this))->BerthelotHercourtEssenRate<CoeffType>::rate_and_derivative(T,rate,drate_dT);
        break;
      case KineticsModel::KOOIJ:
        return (static_cast<const KooijRate<CoeffType>*>(this))->KooijRate<CoeffType>::rate_and_derivative(T,rate,drate_dT);
        break;
      case KineticsModel::VANTHOFF:
        return (static_cast<const VantHoffRate<CoeffType>*>(this))->VantHoffRate<CoeffType>::rate_and_derivative(T,rate,drate_dT);
        break;
      default:
        antioch_error();
        break;
      }
  }

} // end namespace Antioch

#endif
