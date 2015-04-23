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
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef ANTIOCH_REACTION_ENUM_H
#define ANTIOCH_REACTION_ENUM_H

namespace Antioch
{
  namespace ReactionType
  {

    enum ReactionType { ELEMENTARY = 0,
                        DUPLICATE,
                        THREE_BODY, 
                        LINDEMANN_FALLOFF,
                        TROE_FALLOFF,
                        LINDEMANN_FALLOFF_THREE_BODY,
                        TROE_FALLOFF_THREE_BODY};

  enum Parameters{ NOT_FOUND = 0,
                   EFFICIENCIES,
                   TROE_ALPHA,
                   TROE_T1,
                   TROE_T2,
                   TROE_T3
                  };

  } // end namespace ReactionType


} // end namespace Antioch

#endif // ANTIOCH_REACTION_ENUM_H
