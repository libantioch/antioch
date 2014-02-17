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
#include "antioch/antioch_asserts.h"
#include "antioch/kinetics_enum.h"
#include "antioch/metaprogramming.h"

//C++
#include <string>
#include <vector>
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

  template <typename CoeffType, typename VectorCoeffType>
  class PhotochemicalRate;

  /*!\class KineticsType
   * \brief base class for kinetics models
   *
   * Kinetics models are:
   *   - Hercourt Essen (HercourtEssenRate)
   *   - Berthelot (BerthelotRate)
   *   - Arrhenius (ArrheniusRate)
   *   - Berthelot Hercourt Essen (BerthelotHercourtEssenRate)
   *   - Kooij (KooijRate)
   *   - Van't Hoff (VantHoffRate)
   *   - photochemical model (PhotochemicalRate)
   *
   */
  template <typename CoeffType, typename VectorCoeffType = std::vector<CoeffType> >
  class KineticsType{
  public:
    KineticsType(const KineticsModel::KineticsModel type);
    virtual ~KineticsType();

    //!
    template <typename StateType>
    StateType operator()(const StateType& T) const;

    //!
    template <typename StateType>
    StateType derivative( const StateType& T ) const;

    //!
    template <typename StateType>
    void compute_rate_and_derivative(const StateType& T, StateType& rate, StateType& drate_dT) const;

    //!Particles flux supports
    void calculate_rate_constant(const VectorCoeffType &abs,const VectorCoeffType &flux, bool x_update);

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
  template <typename CoeffType, typename VectorCoeffType>
  inline
  void KineticsType<CoeffType,VectorCoeffType>::print(std::ostream& os) const
  {
    os << numeric();
  }

  template <typename CoeffType, typename VectorCoeffType>
  inline
  KineticsType<CoeffType,VectorCoeffType>::KineticsType(const KineticsModel::KineticsModel type):
    my_type(type)
  {
    return;
  }

  template <typename CoeffType, typename VectorCoeffType>
  inline
  KineticsType<CoeffType,VectorCoeffType>::~KineticsType()
  {
    return;
  }


  template <typename CoeffType, typename VectorCoeffType>
  inline
  void KineticsType<CoeffType,VectorCoeffType>::calculate_rate_constant(const VectorCoeffType &abs,const VectorCoeffType &flux, bool x_update)
  {
    switch(my_type) 
      {
      case(KineticsModel::PHOTOCHEM):
        {
          (static_cast<PhotochemicalRate<CoeffType,VectorCoeffType>*>(this))->calculate_rate_constant(abs,flux,x_update);
        }
      break;
      default:
        {
          antioch_error();
        }
      }
  }

  template <typename CoeffType, typename VectorCoeffType>
  template <typename StateType>
  inline
  StateType KineticsType<CoeffType,VectorCoeffType>::operator()(const StateType& T) const
  {
    switch(my_type) 
      {
      case(KineticsModel::HERCOURT_ESSEN):
        {
          return (static_cast<const HercourtEssenRate<CoeffType>*>(this))->rate(T);
        }
        break;

      case(KineticsModel::BERTHELOT):
        {
          return (static_cast<const BerthelotRate<CoeffType>*>(this))->rate(T);
        }
        break;

      case(KineticsModel::ARRHENIUS):
        {
          return (static_cast<const ArrheniusRate<CoeffType>*>(this))->rate(T);
        }
        break;

      case(KineticsModel::BHE):
        {
          return (static_cast<const BerthelotHercourtEssenRate<CoeffType>*>(this))->rate(T);
        }
        break;

      case(KineticsModel::KOOIJ):
        {
          return (static_cast<const KooijRate<CoeffType>*>(this))->rate(T);
        }
        break;

      case(KineticsModel::VANTHOFF):
        {
          return (static_cast<const VantHoffRate<CoeffType>*>(this))->rate(T);
        }
        break;

      case(KineticsModel::PHOTOCHEM):
        {
          return (static_cast<const PhotochemicalRate<CoeffType,VectorCoeffType>*>(this))->rate(T);
        }
        break;

      default:
        {
          antioch_error();
        }

      } // switch(my_type)

    // Dummy
    return zero_clone(T);
  }

  template <typename CoeffType, typename VectorCoeffType>
  template <typename StateType>
  inline
  StateType KineticsType<CoeffType,VectorCoeffType>::derivative( const StateType& T ) const
  {
    switch(my_type) 
      {
      case(KineticsModel::HERCOURT_ESSEN):
        {
          return (static_cast<const HercourtEssenRate<CoeffType>*>(this))->derivative(T);
        }
        break;

      case(KineticsModel::BERTHELOT):
        {
          return (static_cast<const BerthelotRate<CoeffType>*>(this))->derivative(T);
        }
        break;

      case(KineticsModel::ARRHENIUS):
        {
          return (static_cast<const ArrheniusRate<CoeffType>*>(this))->derivative(T);
        }
        break;

      case(KineticsModel::BHE):
        {
          return (static_cast<const BerthelotHercourtEssenRate<CoeffType>*>(this))->derivative(T);
        }
        break;

      case(KineticsModel::KOOIJ):
        {
          return (static_cast<const KooijRate<CoeffType>*>(this))->derivative(T);
        }
        break;

      case(KineticsModel::VANTHOFF):
        {
          return (static_cast<const VantHoffRate<CoeffType>*>(this))->derivative(T);
        }
        break;

      case(KineticsModel::PHOTOCHEM):
        {
          return (static_cast<const PhotochemicalRate<CoeffType,VectorCoeffType>*>(this))->derivative(T);
        }
        break;

      default:
        {
          antioch_error();
        }

      } // switch(my_type)

    // Dummy
    return zero_clone(T);
  }

  template <typename CoeffType, typename VectorCoeffType>
  template <typename StateType>
  inline
  void KineticsType<CoeffType,VectorCoeffType>::compute_rate_and_derivative(const StateType& T, StateType& rate, StateType& drate_dT) const
  {
    switch (my_type) 
      {
      case(KineticsModel::HERCOURT_ESSEN):
        {
          (static_cast<const HercourtEssenRate<CoeffType>*>(this))->rate_and_derivative(T,rate,drate_dT);
        }
        break;

      case(KineticsModel::BERTHELOT):
        {
          (static_cast<const BerthelotRate<CoeffType>*>(this))->rate_and_derivative(T,rate,drate_dT);
        }
        break;

      case(KineticsModel::ARRHENIUS):
        {
          (static_cast<const ArrheniusRate<CoeffType>*>(this))->rate_and_derivative(T,rate,drate_dT);
        }
        break;

      case(KineticsModel::BHE):
        {
          (static_cast<const BerthelotHercourtEssenRate<CoeffType>*>(this))->rate_and_derivative(T,rate,drate_dT);
        }
        break;

      case(KineticsModel::KOOIJ):
        {
          (static_cast<const KooijRate<CoeffType>*>(this))->rate_and_derivative(T,rate,drate_dT);
        }
        break;

      case(KineticsModel::VANTHOFF):
        {
          (static_cast<const VantHoffRate<CoeffType>*>(this))->rate_and_derivative(T,rate,drate_dT);
        }
        break;

      case(KineticsModel::PHOTOCHEM):
        {
          (static_cast<const PhotochemicalRate<CoeffType,VectorCoeffType>*>(this))->rate_and_derivative(T,rate,drate_dT);
        }
        break;

      default:
        {
          antioch_error();
        }

      } // switch(my_type)
    
    return;
  }

} // end namespace Antioch

#endif
