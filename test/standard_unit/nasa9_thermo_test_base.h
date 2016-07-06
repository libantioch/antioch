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

#ifndef ANTIOCH_NASA9_THERMO_TEST_BASE_H
#define ANTIOCH_NASA9_THERMO_TEST_BASE_H

// Antioch
#include "antioch/cmath_shims.h"

namespace AntiochTesting
{
  template<typename Scalar>
  class NASA9ThermoTestBase
  {
  public:

    void init()
    {
      this->init_N2_coeffs_0_200();
      this->init_N2_coeffs_200_1000();
      this->init_N2_coeffs_1000_6000();
      this->init_N2_coeffs_6000_20000();
      this->init_N2_coeffs_20000_99999();
      this->init_all_N2_coeffs();

      this->init_NO2_coeffs_0_200();
      this->init_NO2_coeffs_200_1000();
      this->init_NO2_coeffs_1000_6000();
      this->init_NO2_coeffs_6000_99999();
      this->init_all_NO2_coeffs();
    }

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

    void init_N2_coeffs_0_200()
    {
      _N2_coeffs_0_200.resize(9,0.0L);
      _N2_coeffs_0_200[2] =  3.50080109e+00L;
      _N2_coeffs_0_200[7] = -1.04381074e+03L;
      _N2_coeffs_0_200[8] =  3.09890983e+00L;
    }

    void init_N2_coeffs_200_1000()
    {
      _N2_coeffs_200_1000.resize(9);
      _N2_coeffs_200_1000[0] =  2.21037122e+04L;
      _N2_coeffs_200_1000[1] = -3.81846145e+02L;
      _N2_coeffs_200_1000[2] =  6.08273815e+00L;
      _N2_coeffs_200_1000[3] = -8.53091381e-03L;
      _N2_coeffs_200_1000[4] =  1.38464610e-05L;
      _N2_coeffs_200_1000[5] = -9.62579293e-09L;
      _N2_coeffs_200_1000[6] =  2.51970560e-12L;
      _N2_coeffs_200_1000[7] =  7.10845911e+02L;
      _N2_coeffs_200_1000[8] = -1.07600320e+01L;
    }

    void init_N2_coeffs_1000_6000()
    {
      _N2_coeffs_1000_6000.resize(9);
      _N2_coeffs_1000_6000[0] =  5.87709908e+05L;
      _N2_coeffs_1000_6000[1] = -2.23924255e+03L;
      _N2_coeffs_1000_6000[2] =  6.06694267e+00L;
      _N2_coeffs_1000_6000[3] = -6.13965296e-04L;
      _N2_coeffs_1000_6000[4] =  1.49179819e-07L;
      _N2_coeffs_1000_6000[5] = -1.92309442e-11L;
      _N2_coeffs_1000_6000[6] =  1.06194871e-15L;
      _N2_coeffs_1000_6000[7] =  1.28320618e+04L;
      _N2_coeffs_1000_6000[8] = -1.58663484e+01L;
    }

    void init_N2_coeffs_6000_20000()
    {
      _N2_coeffs_6000_20000.resize(9);
      _N2_coeffs_6000_20000[0] =  8.30971200e+08L;
      _N2_coeffs_6000_20000[1] = -6.42048187e+05L;
      _N2_coeffs_6000_20000[2] =  2.02020507e+02L;
      _N2_coeffs_6000_20000[3] = -3.06501961e-02L;
      _N2_coeffs_6000_20000[4] =  2.48685558e-06L;
      _N2_coeffs_6000_20000[5] = -9.70579208e-11L;
      _N2_coeffs_6000_20000[6] =  1.43751673e-15L;
      _N2_coeffs_6000_20000[7] =  4.93850663e+06L;
      _N2_coeffs_6000_20000[8] = -1.67204791e+03L;
    }

    void init_N2_coeffs_20000_99999()
    {
      _N2_coeffs_20000_99999.resize(9,0.0L);
      _N2_coeffs_20000_99999[2] =  7.27314605e+00L;
      _N2_coeffs_20000_99999[7] = -2.73414039e+04L;
      _N2_coeffs_20000_99999[8] = -2.92594085e+01L;
    }

    void init_all_N2_coeffs()
    {
      _all_standard_N2_coeffs.insert(  _all_standard_N2_coeffs.end(),
                              _N2_coeffs_200_1000.begin(),
                              _N2_coeffs_200_1000.end() );

      _all_standard_N2_coeffs.insert(  _all_standard_N2_coeffs.end(),
                              _N2_coeffs_1000_6000.begin(),
                              _N2_coeffs_1000_6000.end() );

      _all_standard_N2_coeffs.insert(  _all_standard_N2_coeffs.end(),
                              _N2_coeffs_6000_20000.begin(),
                              _N2_coeffs_6000_20000.end() );
    }

    void init_NO2_coeffs_0_200()
    {
      _NO2_coeffs_0_200.resize(9,0.0L);
      _NO2_coeffs_0_200[2] = 4.14754156e+00L;
      _NO2_coeffs_0_200[7] = 2.86078945e+03L;
      _NO2_coeffs_0_200[8] = 5.19771055e+00L;
    }

    void init_NO2_coeffs_200_1000()
    {
      _NO2_coeffs_200_1000.resize(9);
      _NO2_coeffs_200_1000[0] = -5.64204584e+04L;
      _NO2_coeffs_200_1000[1] =  9.63309734e+02L;
      _NO2_coeffs_200_1000[2] = -2.43451851e+00L;
      _NO2_coeffs_200_1000[3] =  1.92776414e-02L;
      _NO2_coeffs_200_1000[4] = -1.87456362e-05L;
      _NO2_coeffs_200_1000[5] =  9.14553637e-09L;
      _NO2_coeffs_200_1000[6] = -1.77766146e-12L;
      _NO2_coeffs_200_1000[7] = -1.54793043e+03L;
      _NO2_coeffs_200_1000[8] =  4.06785541e+01L;
    }

    void init_NO2_coeffs_1000_6000()
    {
      _NO2_coeffs_1000_6000.resize(9);
      _NO2_coeffs_1000_6000[0] =  7.21271176e+05L;
      _NO2_coeffs_1000_6000[1] = -3.83253763e+03L;
      _NO2_coeffs_1000_6000[2] =  1.11395534e+01L;
      _NO2_coeffs_1000_6000[3] = -2.23801440e-03L;
      _NO2_coeffs_1000_6000[4] =  6.54762164e-07L;
      _NO2_coeffs_1000_6000[5] = -7.61120803e-11L;
      _NO2_coeffs_1000_6000[6] =  3.32829926e-15L;
      _NO2_coeffs_1000_6000[7] =  2.50244718e+04L;
      _NO2_coeffs_1000_6000[8] = -4.30507224e+01L;
    }

    void init_NO2_coeffs_6000_99999()
    {
      _NO2_coeffs_6000_99999.resize(9,0.0L);
      _NO2_coeffs_6000_99999[2] =  1.51725546e+02L;
      _NO2_coeffs_6000_99999[7] = -2.44062907e+06L;
      _NO2_coeffs_6000_99999[8] = -1.41879306e+03L;
    }

    void init_all_NO2_coeffs()
    {
      _all_standard_NO2_coeffs.insert(  _all_standard_NO2_coeffs.end(),
                               _NO2_coeffs_200_1000.begin(),
                               _NO2_coeffs_200_1000.end() );

      _all_standard_NO2_coeffs.insert(  _all_standard_NO2_coeffs.end(),
                               _NO2_coeffs_1000_6000.begin(),
                               _NO2_coeffs_1000_6000.end() );
    }

  protected:

    std::vector<Scalar> _N2_coeffs_0_200;
    std::vector<Scalar> _N2_coeffs_200_1000;
    std::vector<Scalar> _N2_coeffs_1000_6000;
    std::vector<Scalar> _N2_coeffs_6000_20000;
    std::vector<Scalar> _N2_coeffs_20000_99999;

    std::vector<Scalar> _NO2_coeffs_0_200;
    std::vector<Scalar> _NO2_coeffs_200_1000;
    std::vector<Scalar> _NO2_coeffs_1000_6000;
    std::vector<Scalar> _NO2_coeffs_6000_99999;


    std::vector<Scalar> _all_standard_N2_coeffs;
    std::vector<Scalar> _all_standard_NO2_coeffs;

  };

} // end namespace AntiochTesting

#endif // ANTIOCH_NASA9_THERMO_TEST_BASE_H
