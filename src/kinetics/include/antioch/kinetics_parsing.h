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

#ifndef ANTIOCH_KINETICS_PARSING_H
#define ANTIOCH_KINETICS_PARSING_H

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/kinetics_type.h"
#include "antioch/hercourtessen_rate.h"
#include "antioch/berthelot_rate.h"
#include "antioch/arrhenius_rate.h"
#include "antioch/berthelothercourtessen_rate.h"
#include "antioch/kooij_rate.h"
#include "antioch/vanthoff_rate.h"
#include "antioch/reaction_enum.h"

/*! cross-road to kinetics model
 */
namespace Antioch
{

  template<typename CoeffType>
  KineticsType<CoeffType>* build_rate( const std::vector<CoeffType> &data,
                                       const KineticsModel::KineticsModel &kin );


  //! We take here the parameters as:
  //  - \f$A\f$, \f$\beta\f$, \f$E_a\f$, \f$D\f$, \f$\mathrm{T_{ref}}\f$, \f$\mathrm{scale}\f$.
  //  The \f$\mathrm{scale}\f$ parameter is the factor of \f$E_a\f$ from its unit
  //  to K.
  //
  //  We check the number of parameters then initialize.
  template<typename CoeffType>
  inline
  KineticsType<CoeffType>* build_rate( const std::vector<CoeffType> &data,
                                       const KineticsModel::KineticsModel &kin )
  {
    using std::pow;

    KineticsType<CoeffType>* rate = NULL;

    switch(kin)
      {
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

      case KineticsModel::KOOIJ:
        {
          antioch_assert_equal_to(5,data.size());
          rate = new KooijRate<CoeffType>(data[0],data[1],data[2],data[3],data[4]);//Cf,eta,Ea,Tref,scale
        }
        break;

      case KineticsModel::VANTHOFF:
        {
          antioch_assert_equal_to(6,data.size());
          rate = new VantHoffRate<CoeffType>(data[0],data[1],data[2],data[3],data[4],data[5]);//Cf,eta,Ea,D,Tref,scale
        }
        break;

      default:
        {
          antioch_error();
        }

      } // switch(kin)

    return rate;
  }

} // end namespace Antioch

#endif // ANTIOCH_REACTION_PARSING_H
