//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014-2016 Paul T. Bauman, Benjamin S. Kirk,
//                         Sylvain Plessis, Roy H. Stonger
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

#ifndef ANTIOCH_CEA_CURVE_FIT_H
#define ANTIOCH_CEA_CURVE_FIT_H

#include "antioch/nasa9_curve_fit.h"

namespace Antioch
{
  /*!\class CEACurveFit

    This is a synomym for NASA9CurveFit, the NASA9
    name ensures consistency with the NASA7 objects
    while the CEA name provides backward compatiblity
    and a more `physics-based' name.

    Note that as nothing happens in the destructor (and
    nothing should ever), no need to get virtual in NASA9.
  */
  template<typename CoeffType=double>
  class CEACurveFit: public NASA9CurveFit<CoeffType>
  {
  public:

    CEACurveFit( const std::vector<CoeffType>& coeffs );

    CEACurveFit( const std::vector<CoeffType>& coeffs, const std::vector<CoeffType> & temps );

    ~CEACurveFit(){}
  };

  template<typename CoeffType>
  inline
  CEACurveFit<CoeffType>::CEACurveFit( const std::vector<CoeffType>& coeffs )
    :NASA9CurveFit<CoeffType>(coeffs)
  {}

  template<typename CoeffType>
  inline
  CEACurveFit<CoeffType>::CEACurveFit( const std::vector<CoeffType>& coeffs, const std::vector<CoeffType> & temps )
    :NASA9CurveFit<CoeffType>(coeffs,temps)
  {}

} // end namespace Antioch

#endif // ANTIOCH_CEA_CURVE_FIT_H
