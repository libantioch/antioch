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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

//This class
#include "antioch/reaction.h"
//Antioch
#include "antioch/elementary_reaction.h"
#include "antioch/duplicate_reaction.h"
#include "antioch/threebody_reaction.h"
#include "antioch/falloff_reaction.h"
#include "antioch/lindemann_falloff.h"
#include "antioch/troe_falloff.h"
#include "antioch/antioch_asserts.h"

namespace Antioch
{

  template<typename CoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  StateType Reaction<CoeffType>::compute_forward_rate_coefficient( const VectorStateType& molar_densities,
                                        const StateType& T) const
  {
     switch(_type)
     {
       case(ReactionType::ELEMENTARY):
        return static_cast<ElementaryReaction<CoeffType>*>(this)->ElementaryReaction<CoeffType>::compute_forward_rate_coefficient(molar_densities,T);
        break;
       case(ReactionType::DUPLICATE):
        return static_cast<DuplicateReaction<CoeffType>*>(this)->DuplicateReaction<CoeffType>::compute_forward_rate_coefficient(molar_densities,T);
        break;
       case(ReactionType::THREE_BODY):
        return static_cast<ThreeBodyReaction<CoeffType>*>(this)->ThreeBodyReaction<CoeffType>::compute_forward_rate_coefficient(molar_densities,T);
        break;
       case(ReactionType::LINDEMANN_FALLOFF):
        return static_cast<FalloffReaction<CoeffType,LindemannFalloff<CoeffType> >*>(this)->FalloffReaction<CoeffType,LindemannFalloff<CoeffType> >::compute_forward_rate_coefficient(molar_densities,T);
        break;
       case(ReactionType::TROE_FALLOFF):
        return static_cast<FalloffReaction<CoeffType,TroeFalloff<CoeffType> >*>(this)->FalloffReaction<CoeffType,TroeFalloff<CoeffType> >::compute_forward_rate_coefficient(molar_densities,T);
        break;
       default:
        antioch_error();
        break;
     }
  }
    
  template<typename CoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  void Reaction<CoeffType>::compute_forward_rate_coefficient_and_derivatives( const VectorStateType& molar_densities,
                                                           const StateType& T, 
                                                           StateType& kfwd, 
                                                           StateType& dkfwd_dT,
                                                           VectorStateType& dkfwd_dY) const
  {
     switch(_type)
     {
       case(ReactionType::ELEMENTARY):
        return static_cast<ElementaryReaction<CoeffType>*>(this)->ElementaryReaction<CoeffType>::compute_forward_rate_coefficient_and_derivatives(molar_densities,T,kfwd,dkfwd_dT,dkfwd_dY);
        break;
       case(ReactionType::DUPLICATE):
        return static_cast<DuplicateReaction<CoeffType>*>(this)->DuplicateReaction<CoeffType>::compute_forward_rate_coefficient_and_derivatives(molar_densities,T,kfwd,dkfwd_dT,dkfwd_dY);
        break;
       case(ReactionType::THREE_BODY):
        return static_cast<ThreeBodyReaction<CoeffType>*>(this)->ThreeBodyReaction<CoeffType>::compute_forward_rate_coefficient_and_derivatives(molar_densities,T,kfwd,dkfwd_dT,dkfwd_dY);
        break;
       case(ReactionType::LINDEMANN_FALLOFF):
        return static_cast<FalloffReaction<CoeffType,LindemannFalloff<CoeffType> >*>(this)->FalloffReaction<CoeffType,LindemannFalloff<CoeffType> >::compute_forward_rate_coefficient_and_derivatives(molar_densities,T,kfwd,dkfwd_dT,dkfwd_dY);
        break;
       case(ReactionType::TROE_FALLOFF):
        return static_cast<FalloffReaction<CoeffType,TroeFalloff<CoeffType> >*>(this)->FalloffReaction<CoeffType,TroeFalloff<CoeffType> >::compute_forward_rate_coefficient_and_derivatives(molar_densities,T,kfwd,dkfwd_dT,dkfwd_dY);
        break;
       default:
        antioch_error();
        break;
     }
  }

} // end namespace Antioch
