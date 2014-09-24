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

#ifndef ANTIOCH_CEA_MIXTURE_H
#define ANTIOCH_CEA_MIXTURE_H

// Antioch
#include "antioch/nasa_mixture.h"

// C++

namespace Antioch
{
  
  // this is now deprecated, but
  // kept for backward compatibility
  // this is but a name for
  // a partial specialization
  template<typename CoeffType=double>
  class CEAThermoMixture: public NASAThermoMixture<CoeffType,CEACurveFit<CoeffType> >
  {
  public:

    CEAThermoMixture( const ChemicalMixture<CoeffType>& chem_mixture ):
        NASAThermoMixture(chem_mixture)
        {antioch_deprecated();}
                
    ~CEAThermoMixture();

  };

} // end namespace Antioch

#endif // ANTIOCH_CEA_MIXTURE_H
