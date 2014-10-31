//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014 Paul T. Bauman, Benjamin S. Kirk, Sylvain Plessis,
//                    Roy H. Stonger
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
#include "antioch/kinetics_conditions.h"

//C++
#include <string>
#include <vector>
#include <iostream>

namespace Antioch{

  template <typename CoeffType>
  class ConstantRate;

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
   *   - constant (ConstantRate)
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

    KineticsModel::KineticsModel type() const;

    //!
    template <typename StateType, typename VectorStateType>
    StateType operator()(const KineticsConditions<StateType,VectorStateType> & conditions) const;

    //!
    template <typename StateType, typename VectorStateType>
    StateType derivative( const KineticsConditions<StateType,VectorStateType> & conditions ) const;

    //!
    template <typename StateType, typename VectorStateType>
    void compute_rate_and_derivative(const KineticsConditions<StateType, VectorStateType>& conditions, StateType& rate, StateType& drate_dT) const;

    //!
    void set_index(unsigned int nr);

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
    unsigned int                 my_index;
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
    my_type(type),
    my_index(0)
  {
    return;
  }

  template <typename CoeffType, typename VectorCoeffType>
  inline
  KineticsModel::KineticsModel KineticsType<CoeffType,VectorCoeffType>::type() const
  {
    return my_type;
  }

  template <typename CoeffType, typename VectorCoeffType>
  inline
  KineticsType<CoeffType,VectorCoeffType>::~KineticsType()
  {
    return;
  }

  template <typename CoeffType, typename VectorCoeffType>
  inline
  void KineticsType<CoeffType,VectorCoeffType>::set_index(unsigned int nr)
  {
    my_index = nr;
  }
    
  template <typename CoeffType, typename VectorCoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  StateType KineticsType<CoeffType,VectorCoeffType>::operator()(const KineticsConditions<StateType,VectorStateType>& conditions) const
  {
    switch(my_type) 
      {
      case(KineticsModel::CONSTANT):
        {
          return (static_cast<const ConstantRate<CoeffType>*>(this))->rate(conditions.T());
        }
        break;

      case(KineticsModel::HERCOURT_ESSEN):
        {
          return (static_cast<const HercourtEssenRate<CoeffType>*>(this))->rate(conditions.T());
        }
        break;

      case(KineticsModel::BERTHELOT):
        {
          return (static_cast<const BerthelotRate<CoeffType>*>(this))->rate(conditions.T());
        }
        break;

      case(KineticsModel::ARRHENIUS):
        {
          return (static_cast<const ArrheniusRate<CoeffType>*>(this))->rate(conditions.T());
        }
        break;

      case(KineticsModel::BHE):
        {
          return (static_cast<const BerthelotHercourtEssenRate<CoeffType>*>(this))->rate(conditions.T());
        }
        break;

      case(KineticsModel::KOOIJ):
        {
          return (static_cast<const KooijRate<CoeffType>*>(this))->rate(conditions.T());
        }
        break;

      case(KineticsModel::VANTHOFF):
        {
          return (static_cast<const VantHoffRate<CoeffType>*>(this))->rate(conditions.T());
        }
        break;

      case(KineticsModel::PHOTOCHEM):
        {
          return (static_cast<const PhotochemicalRate<CoeffType,VectorCoeffType>*>(this))->rate(conditions.particle_flux(my_index));
        }
        break;

      default:
        {
          antioch_error();
        }

      } // switch(my_type)

    // Dummy
    return zero_clone(conditions.T());
  }

  template <typename CoeffType, typename VectorCoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  StateType KineticsType<CoeffType,VectorCoeffType>::derivative( const KineticsConditions<StateType, VectorStateType> & conditions ) const
  {
    switch(my_type) 
      {
      case(KineticsModel::CONSTANT):
        {
          return (static_cast<const ConstantRate<CoeffType>*>(this))->derivative(conditions.T());
        }
        break;

      case(KineticsModel::HERCOURT_ESSEN):
        {
          return (static_cast<const HercourtEssenRate<CoeffType>*>(this))->derivative(conditions.T());
        }
        break;

      case(KineticsModel::BERTHELOT):
        {
          return (static_cast<const BerthelotRate<CoeffType>*>(this))->derivative(conditions.T());
        }
        break;

      case(KineticsModel::ARRHENIUS):
        {
          return (static_cast<const ArrheniusRate<CoeffType>*>(this))->derivative(conditions.T());
        }
        break;

      case(KineticsModel::BHE):
        {
          return (static_cast<const BerthelotHercourtEssenRate<CoeffType>*>(this))->derivative(conditions.T());
        }
        break;

      case(KineticsModel::KOOIJ):
        {
          return (static_cast<const KooijRate<CoeffType>*>(this))->derivative(conditions.T());
        }
        break;

      case(KineticsModel::VANTHOFF):
        {
          return (static_cast<const VantHoffRate<CoeffType>*>(this))->derivative(conditions.T());
        }
        break;

      case(KineticsModel::PHOTOCHEM):
        {
          return (static_cast<const PhotochemicalRate<CoeffType,VectorCoeffType>*>(this))->derivative(conditions.particle_flux(my_index)); 
        }
        break;

      default:
        {
          antioch_error();
        }

      } // switch(my_type)

    // Dummy
    return zero_clone(conditions.T());
  }

  template <typename CoeffType, typename VectorCoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  void KineticsType<CoeffType,VectorCoeffType>::compute_rate_and_derivative(const KineticsConditions<StateType,VectorStateType>& conditions, 
                                                                                  StateType& rate, StateType& drate_dT) const
  {
    switch (my_type) 
      {
      case(KineticsModel::CONSTANT):
        {
          (static_cast<const ConstantRate<CoeffType>*>(this))->rate_and_derivative(conditions.T(),rate,drate_dT);
        }
        break;

      case(KineticsModel::HERCOURT_ESSEN):
        {
          (static_cast<const HercourtEssenRate<CoeffType>*>(this))->rate_and_derivative(conditions.T(),rate,drate_dT);
        }
        break;

      case(KineticsModel::BERTHELOT):
        {
          (static_cast<const BerthelotRate<CoeffType>*>(this))->rate_and_derivative(conditions.T(),rate,drate_dT);
        }
        break;

      case(KineticsModel::ARRHENIUS):
        {
          (static_cast<const ArrheniusRate<CoeffType>*>(this))->rate_and_derivative(conditions.T(),rate,drate_dT);
        }
        break;

      case(KineticsModel::BHE):
        {
          (static_cast<const BerthelotHercourtEssenRate<CoeffType>*>(this))->rate_and_derivative(conditions.T(),rate,drate_dT);
        }
        break;

      case(KineticsModel::KOOIJ):
        {
          (static_cast<const KooijRate<CoeffType>*>(this))->rate_and_derivative(conditions.T(),rate,drate_dT);
        }
        break;

      case(KineticsModel::VANTHOFF):
        {
          (static_cast<const VantHoffRate<CoeffType>*>(this))->rate_and_derivative(conditions.T(),rate,drate_dT);
        }
        break;

      case(KineticsModel::PHOTOCHEM):
        {
          (static_cast<const PhotochemicalRate<CoeffType,VectorCoeffType>*>(this))->rate_and_derivative(conditions.particle_flux(my_index),rate,drate_dT);
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
