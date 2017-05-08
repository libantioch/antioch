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

#ifndef ANTIOCH_NASA_TEST_HELPER_H
#define ANTIOCH_NASA_TEST_HELPER_H

namespace AntiochTesting
{
  template<typename Scalar>
  class NASA7ThermoTestHelper
  {
  public:

    Scalar cp_exact( Scalar T, Scalar a0, Scalar a1, Scalar a2, Scalar a3, Scalar a4 )
    {
      return a0 + a1*T + a2*T*T + a3*T*T*T + a4*(T*T*T*T);
    }

    Scalar h_exact( Scalar T, Scalar a0, Scalar a1, Scalar a2, Scalar a3, Scalar a4, Scalar a5 )
    {
      return a0 + a1/2.0L*T + a2/3.0L*T*T + a3/4.0L*T*T*T + a4/5.0L*(T*T*T*T) + a5/T;
    }

    Scalar s_exact( Scalar T, Scalar a0, Scalar a1, Scalar a2, Scalar a3, Scalar a4, Scalar a6 )
    {
      return a0*std::log(T) + a1*T + a2/2.0L*T*T + a3/3.0L*T*T*T + a4/4.0L*(T*T*T*T) + a6;
    }
  };

  template<typename Scalar>
  class NASA9ThermoTestHelper
  {
  public:

    Scalar cp_exact( Scalar T, Scalar a0, Scalar a1, Scalar a2,
                     Scalar a3, Scalar a4, Scalar a5, Scalar a6 )
    {
      return a0/(T*T) + a1/T + a2 + a3*T + a4*(T*T) + a5*(T*T*T) + a6*(T*T*T*T);
    }

    Scalar h_exact( Scalar T, Scalar a0, Scalar a1, Scalar a2,
                    Scalar a3, Scalar a4, Scalar a5, Scalar a6,
                    Scalar a7 )
    {
      return -a0/(T*T) + a1*std::log(T)/T + a2 + a3*T/2.0L + a4*(T*T)/3.0L
        + a5*(T*T*T)/4.0L + a6*(T*T*T*T)/5.0L + a7/T;
    }

    Scalar s_exact( Scalar T, Scalar a0, Scalar a1, Scalar a2,
                    Scalar a3, Scalar a4, Scalar a5, Scalar a6,
                    Scalar a8 )
    {
      return -a0/(2.L*T*T) - a1/T + a2*std::log(T) + a3*T + a4*(T*T)/2.0L
        + a5*(T*T*T)/3.0L + a6*(T*T*T*T)/4.0L + a8;
    }
  };

} // end namespace AntiochTesting

#endif // ANTIOCH_NASA_TEST_HELPER_H
