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

#ifndef ANTIOCH_NASA7_THERMO_TEST_BASE_H
#define ANTIOCH_NASA7_THERMO_TEST_BASE_H

// Antioch
#include "antioch/cmath_shims.h"

namespace AntiochTesting
{
  template<typename Scalar>
  class NASA7ThermoTestBase
  {
  public:

    virtual void init()
    {
      this->init_H2_coeffs_200_1000();
      this->init_H2_coeffs_1000_3500();
      this->init_all_H2_coeffs();

      this->init_N2_coeffs_300_1000();
      this->init_N2_coeffs_1000_5000();
      this->init_all_N2_coeffs();
    }

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

    void init_H2_coeffs_200_1000()
    {
      _H2_coeffs_200_1000.resize(7);
      _H2_coeffs_200_1000[0] =  2.344331120E+00L;
      _H2_coeffs_200_1000[1] =  7.980520750E-03L;
      _H2_coeffs_200_1000[2] = -1.947815100E-05L;
      _H2_coeffs_200_1000[3] =  2.015720940E-08L;
      _H2_coeffs_200_1000[4] = -7.376117610E-12L;
      _H2_coeffs_200_1000[5] = -9.179351730E+02L;
      _H2_coeffs_200_1000[6] =  6.830102380E-01L;
    }

    void init_H2_coeffs_1000_3500()
    {
      _H2_coeffs_1000_3500.resize(7);
      _H2_coeffs_1000_3500[0] =  3.337279200E+00L;
      _H2_coeffs_1000_3500[1] = -4.940247310E-05L;
      _H2_coeffs_1000_3500[2] =  4.994567780E-07L;
      _H2_coeffs_1000_3500[3] = -1.795663940E-10L;
      _H2_coeffs_1000_3500[4] =  2.002553760E-14L;
      _H2_coeffs_1000_3500[5] = -9.501589220E+02L;
      _H2_coeffs_1000_3500[6] = -3.205023310E+00L;
    }

    void init_all_H2_coeffs()
    {
      _all_standard_H2_coeffs.insert(  _all_standard_H2_coeffs.end(),
                                       _H2_coeffs_200_1000.begin(),
                                       _H2_coeffs_200_1000.end() );

      _all_standard_H2_coeffs.insert(  _all_standard_H2_coeffs.end(),
                                       _H2_coeffs_1000_3500.begin(),
                                       _H2_coeffs_1000_3500.end() );
    }

    void init_N2_coeffs_300_1000()
    {
      _N2_coeffs_300_1000.resize(7);
      _N2_coeffs_300_1000[0] =  3.298677000E+00L;
      _N2_coeffs_300_1000[1] =  1.408240400E-03L;
      _N2_coeffs_300_1000[2] = -3.963222000E-06L;
      _N2_coeffs_300_1000[3] =  5.641515000E-09L;
      _N2_coeffs_300_1000[4] = -2.444854000E-12L;
      _N2_coeffs_300_1000[5] = -1.020899900E+03L;
      _N2_coeffs_300_1000[6] =  3.950372000E+00L;
    }

    void init_N2_coeffs_1000_5000()
    {
      _N2_coeffs_1000_5000.resize(7);
      _N2_coeffs_1000_5000[0] =  2.926640000E+00L;
      _N2_coeffs_1000_5000[1] =  1.487976800E-03L;
      _N2_coeffs_1000_5000[2] = -5.684760000E-07L;
      _N2_coeffs_1000_5000[3] =  1.009703800E-10L;
      _N2_coeffs_1000_5000[4] = -6.753351000E-15L;
      _N2_coeffs_1000_5000[5] = -9.227977000E+02L;
      _N2_coeffs_1000_5000[6] =  5.980528000E+00L;
    }

    void init_all_N2_coeffs()
    {
      _all_standard_N2_coeffs.insert(  _all_standard_N2_coeffs.end(),
                                       _N2_coeffs_300_1000.begin(),
                                       _N2_coeffs_300_1000.end() );

      _all_standard_N2_coeffs.insert(  _all_standard_N2_coeffs.end(),
                                       _N2_coeffs_1000_5000.begin(),
                                       _N2_coeffs_1000_5000.end() );
    }

  protected:

    std::vector<Scalar> _H2_coeffs_200_1000;
    std::vector<Scalar> _H2_coeffs_1000_3500;

    std::vector<Scalar> _N2_coeffs_300_1000;
    std::vector<Scalar> _N2_coeffs_1000_5000;


    std::vector<Scalar> _all_standard_H2_coeffs;
    std::vector<Scalar> _all_standard_N2_coeffs;

  };

} // end namespace AntiochTesting

#endif // ANTIOCH_NASA7_THERMO_TEST_BASE_H
