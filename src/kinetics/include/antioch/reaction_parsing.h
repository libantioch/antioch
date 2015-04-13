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
//
// $Id: reaction_enum.h 38785 2013-04-19 18:46:58Z splessis $
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef ANTIOCH_REACTION_PARSING_H
#define ANTIOCH_REACTION_PARSING_H

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/reaction.h"
#include "antioch/elementary_reaction.h"
#include "antioch/duplicate_reaction.h"
#include "antioch/threebody_reaction.h"
#include "antioch/falloff_reaction.h"
#include "antioch/falloff_threebody_reaction.h"

namespace Antioch
{

  template<typename CoeffType>
  Reaction<CoeffType>* build_reaction( const unsigned int n_species, 
                                       const std::string& equation, 
                                       const bool &reversible,
                                       const ReactionType::ReactionType& type , 
                                       const KineticsModel::KineticsModel& kin );


  template <typename CoeffType>
  void reset_parameter_of_chemical_process(Reaction<CoeffType> * reac,
                                           ReactionType::Parameters parameter,
                                           const CoeffType & new_value,
                                           unsigned int species = 0); // if resetting epsilon

// ----------------------

  template<typename CoeffType>
  inline
  Reaction<CoeffType>* build_reaction( const unsigned int n_species, 
                                       const std::string& equation, 
                                       const bool &reversible,
                                       const ReactionType::ReactionType& type , 
                                       const KineticsModel::KineticsModel& kin )
  {
    Reaction<CoeffType>* reaction = NULL;

    switch(type)
      {
      case(ReactionType::ELEMENTARY):
        {
          reaction = new ElementaryReaction<CoeffType>(n_species,equation,reversible,kin);
        }
        break;

      case(ReactionType::DUPLICATE):
        {
          reaction = new DuplicateReaction<CoeffType>(n_species,equation,reversible,kin);
        }
        break;

      case(ReactionType::THREE_BODY):
        {
          reaction = new ThreeBodyReaction<CoeffType>(n_species,equation,reversible,kin);
        }
        break;

      case(ReactionType::LINDEMANN_FALLOFF):
        {
          reaction = new FalloffReaction<CoeffType,LindemannFalloff<CoeffType> >(n_species,equation,reversible,type,kin);
        }
      case(ReactionType::TROE_FALLOFF):
        {
          reaction = new FalloffReaction<CoeffType,TroeFalloff<CoeffType> >(n_species,equation,reversible,type,kin);
        }
        break;

      case(ReactionType::LINDEMANN_FALLOFF_THREE_BODY):
        {
          reaction = new FalloffThreeBodyReaction<CoeffType,LindemannFalloff<CoeffType> >(n_species,equation,reversible,type,kin);
        }
      case(ReactionType::TROE_FALLOFF_THREE_BODY):
        {
          reaction = new FalloffThreeBodyReaction<CoeffType,TroeFalloff<CoeffType> >(n_species,equation,reversible,type,kin);
        }
        break;

      default:
        {
          antioch_error();
        }

      } // switch(type)
    
    // Dummy
    return reaction;
  }


  template <typename CoeffType>
  inline
  void reset_parameter_of_chemical_process(Reaction<CoeffType> * reac,
                                           ReactionType::Parameters parameter,
                                           const CoeffType & new_value,
                                           unsigned int species)
  {
      switch(parameter)
      {
        case ReactionType::EFFICIENCIES:
        {
          reac->set_efficiency(std::string(),species,new_value); // tests inside, if not three-body or species wrong
        }
          break;
        case ReactionType::TROE_ALPHA:
        {
          antioch_assert( (reac->type() == ReactionType::TROE_FALLOFF) || 
                          (reac->type() == ReactionType::TROE_FALLOFF_THREE_BODY) );
          reac->falloff().set_alpha(new_value);
        }
          break;
        case ReactionType::TROE_T1:
        {
          antioch_assert( (reac->type() == ReactionType::TROE_FALLOFF) || 
                          (reac->type() == ReactionType::TROE_FALLOFF_THREE_BODY) );
          reac->falloff().set_T1(new_value);
        }
          break;
        case ReactionType::TROE_T2:
        {
          antioch_assert( (reac->type() == ReactionType::TROE_FALLOFF) || 
                          (reac->type() == ReactionType::TROE_FALLOFF_THREE_BODY) );
          reac->falloff().set_T2(new_value);
        }
          break;
        case ReactionType::TROE_T3:
        {
          antioch_assert( (reac->type() == ReactionType::TROE_FALLOFF) || 
                          (reac->type() == ReactionType::TROE_FALLOFF_THREE_BODY) );
          reac->falloff().set_T3(new_value);
        }
          break;
        default:
        {
          antioch_error();
        }
          break;
      }
  }


} // end namespace Antioch

#endif // ANTIOCH_REACTION_PARSING_H
