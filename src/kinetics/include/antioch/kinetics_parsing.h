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
  KineticsType<CoeffType> * get_rate_ptr(const std::vector<CoeffType> &data,
                                         const KinMod::KinMod &kin,
                                         const CoeffType T0 = 1.);


//! We take here the parameters as:
//  - \f$A\f$ is not the reduced parameters, thus, in corresponding models \f$A_r = A\times \mathrm{T_{\text{ref}}}^\beta\f$
//  - \f$\beta\f$ is the same
//  - \f$E_a\f$ is the reduced parameter (in K)
//  - \f$D\f$ is the same
//
//  We check the number of parameters then initialize.
  template<typename CoeffType>
  inline
  KineticsType<CoeffType> * get_rate_ptr(const std::vector<CoeffType> &data,
                                         const KinMod::KinMod &kin,
                                         const CoeffType T0)
  {
     using std::pow;
     switch(kin)
     {
       case KinMod::HERCOURT_ESSEN:
       {
         antioch_assert_equal_to(2,data.size());
         HercourtEssenRate<CoeffType> * HErate = new HercourtEssenRate<CoeffType>(data[0]*pow(KinMod::Tref/T0,data[1]),data[1]);
         return static_cast<KineticsType<CoeffType>*> (HErate);
         break;
       }
       case KinMod::BERTHELOT:
       {
         antioch_assert_equal_to(2,data.size());
         BerthelotRate<CoeffType> * Brate = new BerthelotRate<CoeffType>(data[0],data[1]);
         return static_cast<KineticsType<CoeffType>*> (Brate);
         break;
       }
       case KinMod::ARRHENIUS:
       {
         antioch_assert_equal_to(2,data.size());
         ArrheniusRate<CoeffType> * Arate = new ArrheniusRate<CoeffType>(data[0],data[1]);
         return static_cast<KineticsType<CoeffType>*> (Arate);
         break;
       }
       case KinMod::BHE:
       {
         antioch_assert_equal_to(3,data.size());
         BerthelotHercourtEssenRate<CoeffType> * BHErate = new BerthelotHercourtEssenRate<CoeffType>(data[0]*pow(KinMod::Tref/T0,data[1]),data[1],data[2]);
         return static_cast<KineticsType<CoeffType>*> (BHErate);
         break;
       }
       case KinMod::KOOIJ:
       {
         antioch_assert_equal_to(3,data.size());
         KooijRate<CoeffType> * Krate = new KooijRate<CoeffType>(data[0]*pow(KinMod::Tref/T0,data[1]),data[1],data[2]);
         return static_cast<KineticsType<CoeffType>*> (Krate);
         break;
       }
       case KinMod::VANTHOFF:
       {
         antioch_assert_equal_to(4,data.size());
         VantHoffRate<CoeffType> * VHrate = new VantHoffRate<CoeffType>(data[0]*pow(KinMod::Tref/T0,data[1]),data[1],data[2],data[3]);
         return static_cast<KineticsType<CoeffType>*> (VHrate);
         break;
       }
       default:
       {
         antioch_error();
         break;
       }
     }
  }

} // end namespace Antioch

#endif // ANTIOCH_REACTION_PARSING_H
