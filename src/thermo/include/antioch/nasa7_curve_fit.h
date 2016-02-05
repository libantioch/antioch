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


#ifndef ANTIOCH_NASA7_CURVE_FIT_H
#define ANTIOCH_NASA7_CURVE_FIT_H

// Antioch
#include "antioch/nasa_curve_fit_base.h"

namespace Antioch
{
  template<typename CoeffType=double>
  class NASA7CurveFit : public NASACurveFitBase<CoeffType>
  {
  public:

    // default: [300,1000] [1000,5000]
    NASA7CurveFit( const std::vector<CoeffType>& coeffs);

    NASA7CurveFit( const std::vector<CoeffType>& coeffs, const std::vector<CoeffType> & temp );
    ~NASA7CurveFit(){};

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
      @returns the value \f$\frac{\partial\left(\frac{g}{\mathrm{R}T}\right)}{\partial T}
                            = \frac{\partial\left(\frac{h}{\mathrm{R}T} - \frac{s}{R}\right)}{\partial T}\f$
        \f[
            \frac{g}{\mathrm{R}T} = - \frac{a_5}{T^2}   - \frac{a_0}{T}
                                    - \frac{a_1}{2}     - \frac{a_2}{3} T
                                    - \frac{a_3}{4} T^2 - \frac{a_4}{5} T^3
        \f]
     */
    template <typename StateType>
    StateType dh_RT_minus_s_R_dT( const TempCache<StateType>& cache) const;

  };


  /* ------------------------- Inline Functions -------------------------*/

  template<typename CoeffType>
  inline
  NASA7CurveFit<CoeffType>::NASA7CurveFit( const std::vector<CoeffType>& coeffs, const std::vector<CoeffType> & temp )
    : NASACurveFitBase<CoeffType>(coeffs,temp)
  {
    this->_n_coeffs = 7;

    this->check_coeff_size();
    this->check_temp_coeff_size_consistency();
  }

  template<typename CoeffType>
  inline
  NASA7CurveFit<CoeffType>::NASA7CurveFit( const std::vector<CoeffType>& coeffs)
    : NASACurveFitBase<CoeffType>(coeffs,std::vector<CoeffType>())
  {
    this->_n_coeffs = 7;

    this->_temp.resize(3);
    this->_temp[0] = 300.L;
    this->_temp[1] = 1000.L;
    this->_temp[2] = 5000.L;

    this->check_coeff_size();
    this->check_temp_coeff_size_consistency();
  }

  template<typename CoeffType>
  template <typename StateType>
  inline
  const StateType NASA7CurveFit<CoeffType>::cp_over_R(const TempCache<StateType>& cache) const
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
  StateType NASA7CurveFit<CoeffType>::h_over_RT( const TempCache<StateType>& cache) const
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
  StateType NASA7CurveFit<CoeffType>::s_over_R( const TempCache<StateType>& cache) const
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
  NASA7CurveFit<CoeffType>::h_RT_minus_s_R( const TempCache<StateType>& cache) const
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
   StateType NASA7CurveFit<CoeffType>::dh_RT_minus_s_R_dT( const TempCache<StateType>& cache) const
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

#endif //ANTIOCH_NASA7_CURVE_FIT_H
