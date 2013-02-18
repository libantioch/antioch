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
#include <cmath>

// This class
#include "antioch/wilke_mixture.h"

// Antioch
#include "antioch/chemical_species.h"
#include "antioch/chemical_mixture.h"

namespace Antioch
{
  template<class NumericType>
  WilkeMixture<NumericType>::WilkeMixture( const ChemicalMixture<NumericType>& chem_mixture )
    : _chem_mixture(chem_mixture),
      _Mr_Ms_to_the_one_fourth(chem_mixture.n_species()),
      _denom(chem_mixture.n_species())
  {

    for( unsigned int r = 0; r < chem_mixture.n_species(); r++ )
      {
	_Mr_Ms_to_the_one_fourth[r].resize(chem_mixture.n_species());
	_denom[r].resize(chem_mixture.n_species());

	for( unsigned int s = 0; s < chem_mixture.n_species(); s++ )
	  {
	    const NumericType Mr = chem_mixture.M(r);
	    const NumericType Ms = chem_mixture.M(s);

	    _Mr_Ms_to_the_one_fourth[r][s] = std::pow( Mr/Ms, 0.25 );
	    _denom[r][s] = std::sqrt(8.0*(1.0+Ms/Mr));
	  }
      }

    return;
  }

  template<class NumericType>
  WilkeMixture<NumericType>::~WilkeMixture()
  {
    return;
  }

  /* ------------------------- Instantiate ------------------------- */
  template class WilkeMixture<double>;

} // end namespace Antioch
