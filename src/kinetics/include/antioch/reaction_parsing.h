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

namespace Antioch
{

  template<typename CoeffType>
  Reaction<CoeffType> * get_reaction_ptr( const unsigned int n_species, 
                                          const std::string &equation, 
                                          const ReactionType::ReactionType &type , 
                                          const KinMod::KinMod &kin );

  template<typename CoeffType>
  inline
  Reaction<CoeffType> * get_reaction_ptr( const unsigned int n_species, 
                                          const std::string &equation, 
                                          const ReactionType::ReactionType &type , 
                                          const KinMod::KinMod &kin )
  {
     switch(type)
     {
       case ReactionType::ELEMENTARY:
       {
        ElementaryReaction<CoeffType> * Ereac = new ElementaryReaction<CoeffType>(n_species,equation,kin);
        return static_cast<Reaction<CoeffType>*> (Ereac);
        break;
       }
       case ReactionType::DUPLICATE:
       {
        DuplicateReaction<CoeffType> * Dreac = new DuplicateReaction<CoeffType>(n_species,equation,kin);
        return static_cast<Reaction<CoeffType>*> (Dreac);
        break;
       }
       case ReactionType::THREE_BODY:
       {
        ThreeBodyReaction<CoeffType> * TBreac = new ThreeBodyReaction<CoeffType>(n_species,equation,kin);
        return static_cast<Reaction<CoeffType>*> (TBreac);
        break;
       }
       case ReactionType::LINDEMANN_FALLOFF:
       {
       case ReactionType::TROE_FALLOFF:
       {
        FalloffReaction<CoeffType> * Freac = new FalloffReaction<CoeffType>(n_species,equation,type,kin);
        return static_cast<Reaction<CoeffType>*> (Freac);
        break;
       }
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
