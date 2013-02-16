//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// This class
#include "antioch/cea_curve_fit.h"

namespace Antioch
{
  template<class NumericType>
  CEACurveFit<NumericType>::CEACurveFit( const std::vector<NumericType>& coeffs )
    : _n_coeffs(10),
      _coefficients(coeffs)
  {
    return;
  }

  template<class NumericType>
  CEACurveFit<NumericType>::~CEACurveFit()
  {
    return;
  }

  template<class NumericType>
  unsigned int CEACurveFit<NumericType>::interval(const NumericType T) const
  {
    unsigned int interval = -1;

    /* CEA thermodynamic intervals are:
       [200-1,000], [1,000-6,000], [6,000-20,000] K */
    /*! \todo This could be generalized */
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

  /* ------------------------- Instantiate ------------------------- */
  template class CEACurveFit<double>;

} // end namespace Antioch
