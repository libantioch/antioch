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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef ANTIOCH_NASA_CURVE_FIT_H
#define ANTIOCH_NASA_CURVE_FIT_H

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/metaprogramming_decl.h" // Antioch::rebind
#include "antioch/temp_cache.h" 

// C++
#include <vector>

namespace Antioch
{
  template<typename CoeffType=double>
  class NASACurveFit
  {
  public:

    // default: [300,1000] [1000,5000]
    NASACurveFit( const std::vector<CoeffType>& coeffs);

    NASACurveFit( const std::vector<CoeffType>& coeffs, const std::vector<CoeffType> & temp );
    ~NASACurveFit();

    //! The number of intervals for this NASA curve fit
    unsigned int n_intervals() const;

    //! The interval the input temperature lies in
    /*!
      @returns which curve fit interval the input temperature 
      lies in.  The NASA thermodynamic intervals are 
      [200-1,000], [1,000-6,000], [6,000-20,000] K
     */
    template <typename StateType>
    typename Antioch::rebind<StateType, unsigned int>::type
    interval(const StateType& T) const;

    /*!
      @returns the value \f$\frac{Cp}{\mathrm{R}}\f$
        \f[
	   \frac{Cp}{\mathrm{R}} = \frac{a_0}{T^2} + \frac{a_1}{T} + a_2 + a_3 T +
	                            a_4  T^2 + a_5 T^3 + a_6 T^4
        \f]
     */
    template <typename StateType>
    const StateType cp_over_R (const TempCache<StateType> & cache) const;

    /*!
      @returns the value \f$\frac{h}{\mathrm{R}T}\f$
        \f[
            \frac{h}{\mathrm{R}T} = a0  + \frac{a_1}{2} T + \frac{a_2}{3} T^2 
                                    + \frac{a_3}{4} T^3 + \frac{a_4}{5} T^4 
                                    + \frac{a_5}{T}
        \f]
     */
    template <typename StateType>
    StateType h_over_RT( const TempCache<StateType>& cache) const;

    /*!
      @returns the value \f$\frac{s}{\mathrm{R}}\f$
        \f[
            \frac{s}{\mathrm{R}} =  a_0 \ln(T) + a_1 T + \frac{a_2}{2} T^2 
                                    + \frac{a_3}{3} T^3 + \frac{a_4}{4} T^4 
                                    + a_6
        \f]
     */
    template <typename StateType>
    StateType s_over_R( const TempCache<StateType>& cache) const;

    /*!
      @returns the value \f$\frac{g}{\mathrm{R}T} = \frac{h}{\mathrm{R}T} - \frac{s}{R}\f$
        \f[
            \frac{g}{\mathrm{R}T} =   \frac{a_5}{T}      - a0 \ln(T) 
                                    + (a_0 - a_6) 
                                    - \frac{a_1}{2} T    - \frac{a_2}{6} T^2
                                    - \frac{a_3}{12} T^3 - \frac{a_4}{20} T^4 
        \f]
     */
    template <typename StateType>
    StateType h_RT_minus_s_R( const TempCache<StateType>& cache) const;

    /*!
      @returns the value \f$\frac{\partial\left(\frac{g}{\mathrm{R}T}\right)}\frac{\partial T} 
                            = \frac{\partial\left(\frac{h}{\mathrm{R}T} - \frac{s}{R}\right)}{\partial T}\f$
        \f[
            \frac{g}{\mathrm{R}T} = - \frac{a_5}{T^2}   - \frac{a_0}{T}
                                    - \frac{a_1}{2}     - \frac{a_2}{3} T
                                    - \frac{a_3}{4} T^2 - \frac{a_4}{5} T^3
        \f]
     */
    template <typename StateType>
    StateType dh_RT_minus_s_R_dT( const TempCache<StateType>& cache) const;


    //! @returns a pointer to the coefficients in the interval specified.
    /*!   
      The NASA-style equilibrium curve fits are defined in terms of
      _n_coeffs coefficients for each range fit.
    */
    const CoeffType* coefficients(const unsigned int interval) const;

  protected:

    //! The number of coefficients in each interval
    unsigned int _n_coeffs;

    //! The coefficient data
    /*!
      The coeffcients are packed in linear ordering. That is,
      a0-a6 for the first interval, a0-a6 for the second interval,
      and so on.
     */
    std::vector<CoeffType> _coefficients;

    //! The temperatures
    /*!
      The temperature defining the intervals
     */
    std::vector<CoeffType> _temp;
  };


  /* ------------------------- Inline Functions -------------------------*/

  template<typename CoeffType>
  inline
  NASACurveFit<CoeffType>::NASACurveFit( const std::vector<CoeffType>& coeffs, const std::vector<CoeffType> & temp )
    : _n_coeffs(7),
      _coefficients(coeffs),
      _temp(temp)
  {
      // consistency checks
    antioch_assert_equal_to(_coefficients.size()%7,0);
    antioch_assert_equal_to(_temp.size(),_coefficients.size()/_n_coeffs + 1);
    return;
  }

  template<typename CoeffType>
  inline
  NASACurveFit<CoeffType>::NASACurveFit( const std::vector<CoeffType>& coeffs)
    : _n_coeffs(7),
      _coefficients(coeffs)
  {
      // consistency checks
    antioch_assert_equal_to(_coefficients.size(),14);

    _temp.resize(3);
    _temp[0] = 300.L;
    _temp[1] = 1000.L;
    _temp[2] = 5000.L;
    return;
  }


  template<typename CoeffType>
  inline
  NASACurveFit<CoeffType>::~NASACurveFit()
  {
    return;
  }


  template<typename CoeffType>
  template<typename StateType>
  inline
  typename Antioch::rebind<StateType, unsigned int>::type
  NASACurveFit<CoeffType>::interval(const StateType& T) const
  {
    typedef typename 
      Antioch::rebind<StateType, unsigned int>::type UIntType;
    UIntType interval;
    Antioch::zero_clone(interval, T);

    for(unsigned int i = 1; i < _temp.size(); ++i)
    {
        interval = Antioch::if_else
                   (T > _temp[i-1] && T < _temp[i],
                       i - 1,
                       interval); 
    }
    return interval;
  }


  template<typename CoeffType>
  inline
  unsigned int NASACurveFit<CoeffType>::n_intervals() const 
  { return _coefficients.size() / _n_coeffs; }


  template<typename CoeffType>
  inline
  const CoeffType* NASACurveFit<CoeffType>::coefficients(const unsigned int interval) const
  {
    antioch_assert_less( interval, this->n_intervals() );
    antioch_assert_less_equal( _n_coeffs*(interval+1), _coefficients.size() );
    
    return &_coefficients[_n_coeffs*interval];
  }

  template<typename CoeffType>
  template <typename StateType>
  inline
  const StateType NASACurveFit<CoeffType>::cp_over_R(const TempCache<StateType>& cache) const
  {
    typedef typename
      Antioch::rebind<StateType, unsigned int>::type UIntType;
    const UIntType interval = this->interval(cache.T);
    const unsigned int begin_interval = Antioch::min(interval);
    const unsigned int end_interval = Antioch::max(interval)+1;
    
    // FIXME - this needs expression templates to be faster...

    StateType returnval = Antioch::zero_clone(cache.T);
    for (unsigned int i=begin_interval; i != end_interval; ++i)
      {
        const CoeffType * const a =
          this->coefficients(i);
	returnval = Antioch::if_else
	  (interval == i,
	   StateType(a[0] + a[1]*cache.T + a[2]*cache.T2 + a[3]*cache.T3 + a[4]*cache.T4),
	   returnval);
      }

    return returnval;
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType NASACurveFit<CoeffType>::h_over_RT( const TempCache<StateType>& cache) const
  {
    typedef typename
      Antioch::rebind<StateType, unsigned int>::type UIntType;
    const UIntType interval = this->interval(cache.T);
    const unsigned int begin_interval = Antioch::min(interval);
    const unsigned int end_interval = Antioch::max(interval)+1;
    
    StateType returnval = Antioch::zero_clone(cache.T);

    for (unsigned int i=begin_interval; i != end_interval; ++i)
      {
         const CoeffType *a = this->coefficients(interval);
    
         /* h/RT = a0     + a1*T/2 + a2*T^2/3 + a3*T^3/4 + a4*T^4/5 + a5/T */
        returnval = Antioch::if_else
        ( interval == i,
           StateType(  a[0] + 
                       a[1]*cache.T/2.0L + a[2]*cache.T2/3.0L + a[3]*cache.T3/4.0L +
                       a[4]*cache.T4/5.0L + a[5]/cache.T),
           returnval);
       }
       return returnval;
   }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType NASACurveFit<CoeffType>::s_over_R( const TempCache<StateType>& cache) const
  {
    typedef typename
      Antioch::rebind<StateType, unsigned int>::type UIntType;
    const UIntType interval = this->interval(cache.T);
    const unsigned int begin_interval = Antioch::min(interval);
    const unsigned int end_interval = Antioch::max(interval)+1;
    
    StateType returnval = Antioch::zero_clone(cache.T);

    for (unsigned int i=begin_interval; i != end_interval; ++i)
      {
         const CoeffType *a = this->coefficients(interval);
    
    /* s/R = a0*lnT + a1*T   + a2*T^2/2 + a3*T^3/3 + a4*T^4/4 + a6 */
        returnval = Antioch::if_else
        ( interval == i,
           StateType(   a[0]*cache.lnT 
                      + a[1]*cache.T + a[2]*cache.T2/2.0L + a[3]*cache.T3/3.0L 
                      + a[4]*cache.T4/4.0L + a[6]),
           returnval);
       }
       return returnval;
   }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType
  NASACurveFit<CoeffType>::h_RT_minus_s_R( const TempCache<StateType>& cache) const
  {
    typedef typename
      Antioch::rebind<StateType, unsigned int>::type UIntType;
    const UIntType interval = this->interval(cache.T);
    const unsigned int begin_interval = Antioch::min(interval);
    const unsigned int end_interval = Antioch::max(interval)+1;
    
    StateType returnval = Antioch::zero_clone(cache.T);

    for (unsigned int i=begin_interval; i != end_interval; ++i)
      {
         const CoeffType *a = this->coefficients(interval);
    
    /* h/RT =  a[0]     + a[1]*T/2. + a[2]*T2/3. + a[3]*T3/4. + a[4]*T4/5. + a[5]/T,
       s/R  =  a[0]*lnT + a[1]*T    + a[2]*T2/2. + a[3]*T3/3. + a[4]*T4/4. + a[6]   */
        returnval = Antioch::if_else
        ( interval == i,
	   StateType(a[5]/cache.T - a[0]*cache.lnT
                     + a[0] - a[6]
		     - a[1]/2.L*cache.T   - a[2]*cache.T2/6.L
                     - a[3]*cache.T3/12.L - a[4]*cache.T4/20.L),
           returnval);
       }
       return returnval;
   }

   template <typename CoeffType>
   template <typename StateType>
   inline
   StateType NASACurveFit<CoeffType>::dh_RT_minus_s_R_dT( const TempCache<StateType>& cache) const
   {
    typedef typename
      Antioch::rebind<StateType, unsigned int>::type UIntType;
    const UIntType interval = this->interval(cache.T);
    const unsigned int begin_interval = Antioch::min(interval);
    const unsigned int end_interval = Antioch::max(interval)+1;
    
    // FIXME - this needs expression templates to be faster...

    StateType returnval = Antioch::zero_clone(cache.T);

    /* h/RT = -a[0]/T2    + a[1]*lnT/T + a[2]     + a[3]*T/2. + a[4]*T2/3. + a[5]*T3/4. + a[6]*T4/5. + a[8]/T,
       s/R  = -a[0]/T2/2. - a[1]/T     + a[2]*lnT + a[3]*T    + a[4]*T2/2. + a[5]*T3/3. + a[6]*T4/4. + a[9]   */
    for (unsigned int i=begin_interval; i != end_interval; ++i)
      {
        const CoeffType * const a =
          this->coefficients(i);
	returnval = Antioch::if_else
	  (interval == i,
	   StateType(- a[5]/cache.T2     - a[0]/cache.T
		     - a[1]/2.L          - a[2]*cache.T/3.L
                     - a[3]*cache.T2/4.L - a[4]*cache.T3/5.L),
	   returnval);
      }

    return returnval;
        
   }

} // end namespace Antioch

#endif //ANTIOCH_NASA_CURVE_FIT_H
