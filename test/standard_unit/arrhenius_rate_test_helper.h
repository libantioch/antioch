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

#ifndef ANTIOCH_ARRHENIUS_RATE_TEST_HELPER_H
#define ANTIOCH_ARRHENIUS_RATE_TEST_HELPER_H

#include "antioch_config.h"

#ifdef ANTIOCH_HAVE_CPPUNIT

// Antioch
#include "antioch/arrhenius_rate.h"

namespace AntiochTesting
{

  template <typename Scalar>
  class ArrheniusRateTestHelper
  {
  protected:
    Scalar _Cf, _Ea, _R;

    void reset_params( Scalar Cf, Scalar Ea )
    {
      _Cf = Cf;
      _Ea = Ea;
    }

    void reset_params( Scalar Cf, Scalar Ea, Scalar R )
    {
      _Cf = Cf;
      _Ea = Ea;
      _R = R;
    }

    Scalar value( Scalar T )
    {
      using std::exp;
      return _Cf*exp(-_Ea/(_R*T));
    }

    Scalar deriv( Scalar T )
    {
      using std::exp;
      return _Ea/(_R*T*T)*_Cf *exp(-_Ea/(_R*T));
    }

  };

} // end namespace AntiochTesting

#endif // ANTIOCH_HAVE_CPPUNIT

#endif // ANTIOCH_ARRHENIUS_RATE_TEST_HELPER_H
