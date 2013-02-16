//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef ANTIOCH_CEA_CURVE_FIT_H
#define ANTIOCH_CEA_CURVE_FIT_H

// C++
#include <vector>

// Antioch
#include "antioch/antioch_asserts.h"

namespace Antioch
{
  template<class NumericType>
  class CEACurveFit
  {
  public:
    
    CEACurveFit( const std::vector<NumericType>& coeffs );
    ~CEACurveFit();

    //! The number of intervals for this CEA curve fit
    unsigned int n_intervals() const;

    //! The interval the input temperature lies in
    /*!
      @returns which curve fit interval the input temperature 
      lies in.  The CEA thermodynamic intervals are 
      [200-1,000], [1,000-6,000], [6,000-20,000] K
     */
    unsigned int interval(const NumericType T) const;

    
    //! @returns a pointer to the coefficients in the interval specified.
    /*!   
      The CEA-style equilibrium curve fits are defined in terms of
      _n_coeffs coefficients for each range fit.
    */
    const NumericType* coefficients(const unsigned int interval) const;

  protected:

    //! The number of coefficients in each interval
    const unsigned int _n_coeffs;

    //! The coefficient data
    /*!
      The coeffcients are packed in linear ordering. That is,
      a0-a9 for the first interval, a0-a9 for the second interval,
      and so on.
     */
    const std::vector<NumericType> _coefficients;
  };

  /* ------------------------- Inline Functions -------------------------*/

  template<class NumericType>
  inline
  unsigned int CEACurveFit<NumericType>::n_intervals() const 
  { return _coefficients.size() / _n_coeffs; }

  template<class NumericType>
  inline
  const NumericType* CEACurveFit<NumericType>::coefficients(const unsigned int interval) const
  {
    antioch_assert_less( interval, this->n_intervals() );
    antioch_assert_less_equal( _n_coeffs*(interval+1), _coefficients.size() );
    
    return &_coefficients[_n_coeffs*interval];
  }

} // end namespace Antioch

#endif //ANTIOCH_CEA_CURVE_FIT_H
