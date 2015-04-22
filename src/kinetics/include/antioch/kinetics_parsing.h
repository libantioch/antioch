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

#ifndef ANTIOCH_KINETICS_PARSING_H
#define ANTIOCH_KINETICS_PARSING_H

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/metaprogramming_decl.h"
#include "antioch/physical_constants.h"
#include "antioch/kinetics_type.h"
#include "antioch/constant_rate.h"
#include "antioch/hercourtessen_rate.h"
#include "antioch/berthelot_rate.h"
#include "antioch/arrhenius_rate.h"
#include "antioch/berthelothercourtessen_rate.h"
#include "antioch/kooij_rate.h"
#include "antioch/vanthoff_rate.h"
#include "antioch/photochemical_rate.h"
#include "antioch/reaction_enum.h"

namespace Antioch
{
/*! cross-road to kinetics model
 */

  template<typename CoeffType, typename VectorCoeffType>
  KineticsType<CoeffType, VectorCoeffType>* build_rate( const VectorCoeffType &data,
                                                        const KineticsModel::KineticsModel &kin );

  template<typename CoeffType, typename VectorCoeffType, typename VectorType>
  void reset_rate( KineticsType<CoeffType,VectorCoeffType> & kin, const VectorType & coefs);

  template <typename CoeffType, typename VectorCoeffType>
  void reset_parameter_of_rate(KineticsType<CoeffType,VectorCoeffType> & rate,
                               KineticsModel::Parameters parameter,
                               const CoeffType new_value, const std::string & unit = "SI");

  // vectorized parameter
  template <typename CoeffType, typename VectorCoeffType>
  void reset_parameter_of_rate(KineticsType<CoeffType,VectorCoeffType> & rate,
                               KineticsModel::Parameters parameter,
                               const CoeffType new_value, int l, const std::string & unit = "SI");


//----------------------------------------

  //! We take here the parameters as:
  //  - \f$A\f$, \f$\beta\f$, \f$E_a\f$, \f$D\f$, \f$\mathrm{T_{ref}}\f$, \f$\mathrm{scale}\f$.
  //  The \f$\mathrm{scale}\f$ parameter is the factor of \f$E_a\f$ from its unit
  //  to K.
  //
  //  We check the number of parameters then initialize.
  template<typename CoeffType, typename VectorCoeffType>
  inline
  KineticsType<CoeffType, VectorCoeffType>* build_rate( const VectorCoeffType &data,
                                                        const KineticsModel::KineticsModel &kin )
  {
    using std::pow;

    KineticsType<CoeffType,VectorCoeffType>* rate = NULL;

    switch(kin)
      {
      case(KineticsModel::CONSTANT):
        {
          antioch_assert_equal_to(1,data.size());
          rate = new ConstantRate<CoeffType>(data[0]);//Cf
        }
        break;

      case(KineticsModel::HERCOURT_ESSEN):
        {
          antioch_assert_equal_to(3,data.size());
          rate = new HercourtEssenRate<CoeffType>(data[0],data[1],data[2]);//Cf,eta,Tref
        }
        break;

      case(KineticsModel::BERTHELOT):
        {
          antioch_assert_equal_to(2,data.size());
          rate = new BerthelotRate<CoeffType>(data[0],data[1]);// Cf, D
        }
        break;

      case(KineticsModel::ARRHENIUS):
        {
          antioch_assert_equal_to(3,data.size());
          rate = new ArrheniusRate<CoeffType>(data[0],data[1],data[2]);//Cf,Ea,scale
        }
        break;

      case(KineticsModel::BHE):
        {
          antioch_assert_equal_to(4,data.size());
          rate = new BerthelotHercourtEssenRate<CoeffType>(data[0],data[1],data[2],data[3]);//Cf,eta,D,Tref
        }
        break;

      case(KineticsModel::KOOIJ):
        {
          antioch_assert_equal_to(5,data.size());
          rate = new KooijRate<CoeffType>(data[0],data[1],data[2],data[3],data[4]);//Cf,eta,Ea,Tref,scale
        }
        break;

      case(KineticsModel::VANTHOFF):
        {
          antioch_assert_equal_to(6,data.size());
          rate = new VantHoffRate<CoeffType>(data[0],data[1],data[2],data[3],data[4],data[5]);//Cf,eta,Ea,D,Tref,scale
        }
        break;

      case(KineticsModel::PHOTOCHEM):
        {
          antioch_assert_equal_to(0,data.size()%2); //two vectors in one: lambda then cross-section
          VectorCoeffType cs;
          VectorCoeffType lambda;
          for(unsigned int i = 0; i < data.size()/2; i++)
          {
             lambda.push_back(data[i]);
             cs.push_back(data[i + data.size()/2]); 
          }
          rate = new PhotochemicalRate<CoeffType,VectorCoeffType>(cs,lambda);//cs,lambda
        }
        break;

      default:
        {
          antioch_error();
        }

      } // switch(kin)

    return rate;
  }

  template<typename CoeffType, typename VectorCoeffType, typename VectorType>
  void reset_rate( KineticsType<CoeffType,VectorCoeffType> & kin, const VectorType & coefs)
  {

    switch(kin.type())
      {
      case(KineticsModel::CONSTANT):
        {
          static_cast<ConstantRate<CoeffType>*>(&kin)->reset_coefs(coefs);
        }
        break;

      case(KineticsModel::HERCOURT_ESSEN):
        {
          static_cast< HercourtEssenRate<CoeffType>*>(&kin)->reset_coefs(coefs);
        }
        break;

      case(KineticsModel::BERTHELOT):
        {
          static_cast< BerthelotRate<CoeffType>*>(&kin)->reset_coefs(coefs);
        }
        break;

      case(KineticsModel::ARRHENIUS):
        {
          static_cast< ArrheniusRate<CoeffType>*>(&kin)->reset_coefs(coefs);
        }
        break;

      case(KineticsModel::BHE):
        {
          static_cast< BerthelotHercourtEssenRate<CoeffType>*>(&kin)->reset_coefs(coefs);
        }
        break;

      case(KineticsModel::KOOIJ):
        {
          static_cast< KooijRate<CoeffType>*>(&kin)->reset_coefs(coefs);
        }
        break;

      case(KineticsModel::VANTHOFF):
        {
          static_cast< VantHoffRate<CoeffType>*>(&kin)->reset_coefs(coefs);
        }
        break;

      case(KineticsModel::PHOTOCHEM):
        {
          static_cast<PhotochemicalRate<CoeffType,VectorCoeffType>*>(&kin)->reset_coefs(coefs);
        }
        break;

      default:
        {
          antioch_error();
        }

      } // switch(kin.type())
  }


  template <typename CoeffType, typename VectorCoeffType>
  void reset_parameter_of_rate(KineticsType<CoeffType,VectorCoeffType> & rate,
                               KineticsModel::Parameters parameter,
                               const CoeffType new_value, const std::string & unit)
  {

// this is crude at the moment, no test
// this will be replaced by a unit manager
// at some point, to be able to have an explicit
// custom internal unit system with appropriate testing
    CoeffType new_coef = (unit == "SI")?new_value:
                                        new_value * Units<typename value_type<CoeffType>::type>(unit).get_SI_factor();

// Ea management, we want K, two possibilities now
// 1 - Ea is already in K
// 2 - Ea is in J.mol-1
   if(parameter == KineticsModel::Parameters::E)
   {
      if(unit != "K")
      {
         new_coef = new_coef / Constants::R_universal<typename value_type<CoeffType>::type>();
      }
   }
   

    switch(rate.type())
      {
      case(KineticsModel::CONSTANT):
        {
          static_cast<ConstantRate<CoeffType>*>(&rate)->set_parameter(parameter,new_coef);
        }
        break;

      case(KineticsModel::HERCOURT_ESSEN):
        {
          static_cast< HercourtEssenRate<CoeffType>*>(&rate)->set_parameter(parameter,new_coef);
        }
        break;

      case(KineticsModel::BERTHELOT):
        {
          static_cast< BerthelotRate<CoeffType>*>(&rate)->set_parameter(parameter,new_coef);
        }
        break;

      case(KineticsModel::ARRHENIUS):
        {
          static_cast< ArrheniusRate<CoeffType>*>(&rate)->set_parameter(parameter,new_coef);
        }
        break;

      case(KineticsModel::BHE):
        {
          static_cast< BerthelotHercourtEssenRate<CoeffType>*>(&rate)->set_parameter(parameter,new_coef);
        }
        break;

      case(KineticsModel::KOOIJ):
        {
          static_cast< KooijRate<CoeffType>*>(&rate)->set_parameter(parameter,new_coef);
        }
        break;

      case(KineticsModel::VANTHOFF):
        {
          static_cast< VantHoffRate<CoeffType>*>(&rate)->set_parameter(parameter,new_coef);
        }
        break;

      default:
        {
          antioch_error();
        }

      } // switch(kin.type())
  }

  template <typename CoeffType, typename VectorCoeffType>
  void reset_parameter_of_rate(KineticsType<CoeffType,VectorCoeffType> & rate,
                               KineticsModel::Parameters parameter,
                               const CoeffType new_value, int l, const std::string & unit)
  {

// this is crude at the moment, no test
// this will be replaced by a unit manager
// at some point, to be able to have an explicit
// custom internal unit system with appropriate testing
    CoeffType new_coef = (unit == "SI")?new_value:
                                        new_value * Units<typename value_type<CoeffType>::type>(unit).get_SI_factor();

    switch(rate.type())
    {
      case(KineticsModel::PHOTOCHEM):
        {
          static_cast<PhotochemicalRate<CoeffType,VectorCoeffType>*>(&rate)->set_parameter(parameter,l,new_coef);
        }
        break;

      default:
        {
          antioch_error();
        }

      } // switch(kin.type())
  }



} // end namespace Antioch

#endif // ANTIOCH_REACTION_PARSING_H
