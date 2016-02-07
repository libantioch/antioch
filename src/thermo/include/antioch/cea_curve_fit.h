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

    This class only differs from NASA9CurveFit in the construction. Here,
    we assume that there are 10 coefficients, with the 7th being zero. This is
    exactly the format output from NASA's CEA program and, hence, this class
    was build to enable this compatiblity. Internally, the coefficients are
    remapped to 9 coefficients and is then functionally identical to NASA9CurveFit.
  */
  template<typename CoeffType=double>
  class CEACurveFit: public NASA9CurveFit<CoeffType>
  {
  public:

    CEACurveFit( const std::vector<CoeffType>& coeffs );

    CEACurveFit( const std::vector<CoeffType>& coeffs, const std::vector<CoeffType> & temps );

    ~CEACurveFit(){}

  private:

    void remap_coeffs( const std::vector<CoeffType>& coeffs );

  };

  template<typename CoeffType>
  inline
  CEACurveFit<CoeffType>::CEACurveFit( const std::vector<CoeffType>& coeffs )
    :NASA9CurveFit<CoeffType>()
  {

    if( this->_coefficients.size()%10 != 0 )
      antioch_error_msg("ERROR: Expected CEA style of input for coefficients! Must be a multiple of 10!");

    this->_n_coeffs = 9;

    // If no temp is provided, we assume the standard CEA form.
    this->init_nasa9_temps( coeffs, 10 );

    this->remap_coeffs(coeffs);

    this->check_coeff_size();
    this->check_temp_coeff_size_consistency();
  }

  template<typename CoeffType>
  inline
  CEACurveFit<CoeffType>::CEACurveFit( const std::vector<CoeffType>& coeffs,
                                       const std::vector<CoeffType> & temp )
    :NASA9CurveFit<CoeffType>()
  {
    if( this->_coefficients.size()%10 != 0 )
      antioch_error_msg("ERROR: Expected CEA style of input for coefficients! Must be a multiple of 10!");

    this->_n_coeffs = 9;

    this->_temp = temp;

    this->remap_coeffs(coeffs);

    this->check_coeff_size();
    this->check_temp_coeff_size_consistency();
  }

  template<typename CoeffType>
  inline
  void CEACurveFit<CoeffType>::remap_coeffs( const std::vector<CoeffType>& coeffs )
  {
    this->_coefficients.resize(this->_n_coeffs*(this->_temp.size()-1),0.0);

    for( unsigned int t = 0; t < this->_temp.size()-1; t++ )
      {
        for( unsigned int c = 0; c < 7; c++ )
          {
            unsigned int i = 10*t + c;
            unsigned int j = 9*t + c;
            this->_coefficients[j] = coeffs[i];
          }

        for( unsigned int c = 7; c < 9; c++ )
          {
            unsigned int i = 10*t + c;
            unsigned int j = 9*t + c;
            this->_coefficients[j] = coeffs[i+1];
          }
      }
  }

} // end namespace Antioch

#endif // ANTIOCH_CEA_CURVE_FIT_H
