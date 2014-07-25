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

#ifndef ANTIOCH_MATH_CONSTANTS_H
#define ANTIOCH_MATH_CONSTANTS_H

#include <cmath>

namespace Antioch
{
  namespace Constants
  {
    /*!
     * 1/ln(10)
     */
    template<typename CoeffType>
    inline
    CoeffType log10_to_log()
    {
      return 1.0L/std::log(10.0L);
    }

    /*! pi
     *
     */
    template<typename CoeffType>
    inline
    CoeffType pi()
    {
      static const CoeffType my_pi = std::atan(CoeffType(1))*4;
      return my_pi;
    }

  } // end namespace Constants
} // end namespace Antioch

#endif //ANTIOCH_MATH_CONSTANTS_H
