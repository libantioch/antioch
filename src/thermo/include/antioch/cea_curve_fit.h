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

#ifndef ANTIOCH_CEA_CURVE_FIT_H
#define ANTIOCH_CEA_CURVE_FIT_H

// Antioch
#include "antioch/antioch_asserts.h"

// C++
#include <vector>

namespace Antioch
{
  template<typename CoeffType=double>
  class CEACurveFit
  {
  public:
    
    CEACurveFit( const std::vector<CoeffType>& coeffs );
    ~CEACurveFit();

    //! The number of intervals for this CEA curve fit
    unsigned int n_intervals() const;

    //! The interval the input temperature lies in
    /*!
      @returns which curve fit interval the input temperature 
      lies in.  The CEA thermodynamic intervals are 
      [200-1,000], [1,000-6,000], [6,000-20,000] K
     */
    template <typename StateType>
    unsigned int interval(const StateType& T) const;

    
    //! @returns a pointer to the coefficients in the interval specified.
    /*!   
      The CEA-style equilibrium curve fits are defined in terms of
      _n_coeffs coefficients for each range fit.
    */
    const CoeffType* coefficients(const unsigned int interval) const;

  protected:

    //! The number of coefficients in each interval
    const unsigned int _n_coeffs;

    //! The coefficient data
    /*!
      The coeffcients are packed in linear ordering. That is,
      a0-a9 for the first interval, a0-a9 for the second interval,
      and so on.
     */
    const std::vector<CoeffType> _coefficients;
  };


  /* ------------------------- Inline Functions -------------------------*/

  template<typename CoeffType>
  inline
  CEACurveFit<CoeffType>::CEACurveFit( const std::vector<CoeffType>& coeffs )
    : _n_coeffs(10),
      _coefficients(coeffs)
  {
    return;
  }


  template<typename CoeffType>
  inline
  CEACurveFit<CoeffType>::~CEACurveFit()
  {
    return;
  }


  template<typename CoeffType>
  template<typename StateType>
  inline
  unsigned int CEACurveFit<CoeffType>::interval(const StateType& T) const
  {
    unsigned int interval = -1;

    /* CEA thermodynamic intervals are:
       [200-1,000], [1,000-6,000], [6,000-20,000] K */
    /*! \todo This could be generalized */
    /*! \todo This needs to be vectorizable */
    if (T > 6000.)	  
      {
	interval = 2;
      }
    else if (T > 1000.)
      {
	interval =  1;
      }
    else
      {
	interval = 0;
      }

    return interval;
  }


  template<typename CoeffType>
  inline
  unsigned int CEACurveFit<CoeffType>::n_intervals() const 
  { return _coefficients.size() / _n_coeffs; }


  template<typename CoeffType>
  inline
  const CoeffType* CEACurveFit<CoeffType>::coefficients(const unsigned int interval) const
  {
    antioch_assert_less( interval, this->n_intervals() );
    antioch_assert_less_equal( _n_coeffs*(interval+1), _coefficients.size() );
    
    return &_coefficients[_n_coeffs*interval];
  }

} // end namespace Antioch

#endif //ANTIOCH_CEA_CURVE_FIT_H
