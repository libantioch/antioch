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

// This class
//#include "antioch/kinetics_type.h"
// Antioch
/*#include "antioch/hercourtessen_rate.h"
#include "antioch/berthelot_rate.h"
#include "antioch/arrhenius_rate.h"
#include "antioch/berthelothercourtessen_rate.h"
#include "antioch/kooij_rate.h"
#include "antioch/vanthoff_rate.h"
#include "antioch/antioch_asserts.h"
*/
namespace Antioch
{
/*
  template <typename CoeffType>
  template <typename StateType>
  inline
  StateType KineticsType<CoeffType>::operator()(const StateType& T) const
  {
      switch (my_type) 
      {
      case KinMod::HERCOURT_ESSEN:
        return (static_cast<HercourtEssenRate<CoeffType>*>(this))->HercourtEssenRate<CoeffType>::rate(T);
        break;
      case KinMod::BERTHELOT:
        return (static_cast<BerthelotRate<CoeffType>*>(this))->BerthelotRate<CoeffType>::rate(T);
        break;
      case KinMod::ARRHENIUS:
        return (static_cast<ArrheniusRate<CoeffType>*>(this))->ArrheniusRate<CoeffType>::rate(T);
        break;
      case KinMod::BHE:
        return (static_cast<BerthelotHercourtEssenRate<CoeffType>*>(this))->BerthelotHercourtEssenRate<CoeffType>::rate(T);
        break;
      case KinMod::KOOIJ:
        return (static_cast<KooijRate<CoeffType>*>(this))->KooijRate<CoeffType>::rate(T);
        break;
      case KinMod::VANTHOFF:
        return (static_cast<VantHoffRate<CoeffType>*>(this))->VantHoffRate<CoeffType>::rate(T);
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
      case KinMod::HERCOURT_ESSEN:
        return (static_cast<HercourtEssenRate<CoeffType>*>(this))->HercourtEssenRate<CoeffType>::derivative(T);
        break;
      case KinMod::BERTHELOT:
        return (static_cast<BerthelotRate<CoeffType>*>(this))->BerthelotRate<CoeffType>::derivative(T);
        break;
      case KinMod::ARRHENIUS:
        return (static_cast<ArrheniusRate<CoeffType>*>(this))->ArrheniusRate<CoeffType>::derivative(T);
        break;
      case KinMod::BHE:
        return (static_cast<BerthelotHercourtEssenRate<CoeffType>*>(this))->BerthelotHercourtEssenRate<CoeffType>::derivative(T);
        break;
      case KinMod::KOOIJ:
        return (static_cast<KooijRate<CoeffType>*>(this))->KooijRate<CoeffType>::derivative(T);
        break;
      case KinMod::VANTHOFF:
        return (static_cast<VantHoffRate<CoeffType>*>(this))->VantHoffRate<CoeffType>::derivative(T);
        break;
      default:
        antioch_error();
        break;
      }
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void KineticsType<CoeffType>::rate_and_derivative(const StateType& T, StateType& rate, StateType& drate_dT) const
  {
      switch (my_type) 
      {
      case KinMod::HERCOURT_ESSEN:
        return (static_cast<HercourtEssenRate<CoeffType>*>(this))->HercourtEssenRate<CoeffType>::rate_and_derivative(T,rate,drate_dT);
        break;
      case KinMod::BERTHELOT:
        return (static_cast<BerthelotRate<CoeffType>*>(this))->BerthelotRate<CoeffType>::rate_and_derivative(T,rate,drate_dT);
        break;
      case KinMod::ARRHENIUS:
        return (static_cast<ArrheniusRate<CoeffType>*>(this))->ArrheniusRate<CoeffType>::rate_and_derivative(T,rate,drate_dT);
        break;
      case KinMod::BHE:
        return (static_cast<BerthelotHercourtEssenRate<CoeffType>*>(this))->BerthelotHercourtEssenRate<CoeffType>::rate_and_derivative(T,rate,drate_dT);
        break;
      case KinMod::KOOIJ:
        return (static_cast<KooijRate<CoeffType>*>(this))->KooijRate<CoeffType>::rate_and_derivative(T,rate,drate_dT);
        break;
      case KinMod::VANTHOFF:
        return (static_cast<VantHoffRate<CoeffType>*>(this))->VantHoffRate<CoeffType>::rate_and_derivative(T,rate,drate_dT);
        break;
      default:
        antioch_error();
        break;
      }
  }
*/
} // end namespace Antioch
