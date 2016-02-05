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

#ifndef ANTIOCH_NASA_CURVE_FIT_BASE_H
#define ANTIOCH_NASA_CURVE_FIT_BASE_H

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/metaprogramming_decl.h" // Antioch::rebind
#include "antioch/temp_cache.h"

// C++
#include <vector>

namespace Antioch
{
  template<typename CoeffType=double>
  class NASACurveFitBase
  {
  public:

    NASACurveFitBase( const std::vector<CoeffType>& coeffs, const std::vector<CoeffType>& temp );

    ~NASACurveFitBase(){};

    //! The number of intervals for this NASA9 curve fit
    unsigned int n_intervals() const;

    //! The interval the input temperature lies in
    /*!
      @returns which curve fit interval the input temperature
      lies in.
     */
    template <typename StateType>
    typename Antioch::rebind<StateType, unsigned int>::type
    interval(const StateType& T) const;

    //! @returns a pointer to the coefficients in the interval specified.
    /*!
      The ordering/packing of the coefficients will depend on the subclass
    */
    const CoeffType* coefficients(const unsigned int interval) const;

  protected:

    //! The number of coefficients in each interval
    unsigned int _n_coeffs;

    //! The coefficient data
    /*!
      The coeffcients are packed in linear ordering. That is,
      a0-a9 for the first interval, a0-a9 for the second interval,
      and so on.
     */
    std::vector<CoeffType> _coefficients;

    //! The temperatures
    /*!
      The temperature defining the intervals
     */
    std::vector<CoeffType> _temp;

  private:

    NASACurveFitBase();

  };

  template<typename CoeffType>
  inline
  NASACurveFitBase<CoeffType>::NASACurveFitBase( const std::vector<CoeffType>& coeffs,
                                                 const std::vector<CoeffType>& temp )
    : _n_coeffs(0),
      _coefficients(coeffs),
      _temp(temp)
  {}

  template<typename CoeffType>
  inline
  unsigned int NASACurveFitBase<CoeffType>::n_intervals() const
  { return _coefficients.size() / _n_coeffs; }

  template<typename CoeffType>
  template<typename StateType>
  inline
  typename Antioch::rebind<StateType, unsigned int>::type
  NASACurveFitBase<CoeffType>::interval(const StateType& T) const
  {
    typedef typename
      Antioch::rebind<StateType, unsigned int>::type UIntType;
    UIntType interval;
    Antioch::zero_clone(interval, T);

    for(unsigned int i = 1; i < _temp.size(); ++i)
    {
        interval = Antioch::if_else
          (T > _temp[i-1] && T < _temp[i],
           Antioch::constant_clone(interval, i-1),
           interval );
    }

    return interval;
  }

  template<typename CoeffType>
  inline
  const CoeffType* NASACurveFitBase<CoeffType>::coefficients(const unsigned int interval) const
  {
    antioch_assert_less( interval, this->n_intervals() );
    antioch_assert_less_equal( _n_coeffs*(interval+1), _coefficients.size() );

    return &_coefficients[_n_coeffs*interval];
  }

} // end namespace Antioch

#endif // ANTIOCH_NASA_CURVE_FIT_BASE_H
