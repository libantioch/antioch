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

#ifndef ANTIOCH_PHYSICS_PLACEHOLDER_H
#define ANTIOCH_PHYSICS_PLACEHOLDER_H

// Antioch
#include "antioch/physics_metaprogramming_decl.h"

// C++

namespace Antioch
{
  /*
     An empty all default physics class
   */
  class PhysicsPlaceholder
  {
  public:

      PhysicsPlaceholder(){}
      ~PhysicsPlaceholder(){}
  };

// physics tag, all defaults, nothing happens
   template <>
   struct physical_tag< PhysicsPlaceholder >:
        public physical_tag_base<PhysicsPlaceholder >
   {};

} // namespace Antioch
#endif // ANTIOCH_PHYSICS_PLACEHOLDER_H
