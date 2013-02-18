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

// C++
#include <sstream>

// These functions
#include "antioch/blottner_parsing.h"

// Antioch
#include "antioch/input_utils.h"
#include "antioch/antioch_asserts.h"

namespace Antioch
{
  /* ------------------------- Instantiate ------------------------- */
  template void read_blottner_data_ascii<double>( MixtureViscosity<BlottnerViscosity<double>,double >&,
						  std::istream &in );

  template void read_blottner_data_ascii_default<double>( MixtureViscosity<BlottnerViscosity<double>,double>& mu );

} // end namespace Antioch
