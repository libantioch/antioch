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

#ifndef ANTIOCH_WILKE_MIXTURE_H
#define ANTIOCH_WILKE_MIXTURE_H

// C++
#include <vector>

namespace Antioch
{
  // Forward declarations
  template<class NumericType>
  class ChemicalMixture;

  template<class NumericType>
  class WilkeMixture
  {
  public:

    WilkeMixture( const ChemicalMixture<NumericType>& chem_mixture );
    ~WilkeMixture();

  protected:

    const ChemicalMixture<NumericType>& _chem_mixture;

    //! Cache for numerator term
    /*! \todo We should use a more efficient data structure */
    std::vector<std::vector<NumericType> > _Mr_Ms_to_the_one_fourth;
    
    //! Cache for denominator term
    /*! \todo We should use a more efficient data structure */
    std::vector<std::vector<NumericType> > _denom;

  };

} // end namespace Antioch

#endif // ANTIOCH_WILKE_MIXTURE_H
