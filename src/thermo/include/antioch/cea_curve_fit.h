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

#ifndef ANTIOCH_NASA9_CURVE_FIT_H
#define ANTIOCH_NASA9_CURVE_FIT_H

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/metaprogramming_decl.h" // Antioch::rebind
#include "antioch/temp_cache.h" 

// C++
#include <vector>

namespace Antioch
{
  /*! \class NASA9CurveFit
   *
   *  This class stores the CEA polynomial fit to
   *  the thermodynamics quantities \f$\frac{C_p}{\mathrm{R}}\f$
   *  \f$\frac{h}{\mathrm{R}T}\f$ and \f$\frac{s}{\mathrm{R}T}\f$.
   *  This formulation requires nine coefficients, from \f$a_0\f$
   *  to \f$a_9\f$, \f$a_7\f$ is always equals to zero.
   *
   *  The temperature intervals are imposed as: 
   *   [200--1,000], [1,000--6,000], [6,000--20,000] K.
   *
   *  The equations are:
   *    \f[
   *       \frac{Cp}{\mathrm{R}} = \frac{a_0}{T^2} + \frac{a_1}{T} + a_2 + a_3 T +
   *                                a_4  T^2 + a_5 T^3 + a_6 T^4
   *    \f]
   *
   *    \f[
   *        \frac{h}{\mathrm{R}T} = -\frac{a_0}{T^2} + \frac{a_1}{T} \ln(T) 
   *                                + a2  + \frac{a_3}{2} T + \frac{a_4}{3} T^2 
   *                                + \frac{a_5}{4} T^3 + \frac{a_6}{5} T^4 
   *                                + \frac{a_8}{T}
   *    \f]
   *
   *    \f[
   *        \frac{s}{\mathrm{R}} = -\frac{a_0}{2T^2} - \frac{a_1}{T}
   *                                + a_2 \ln(T) + a_3 T + \frac{a_4}{2} T^2 
   *                                + \frac{a_5}{3} T^3 + \frac{a_6}{4} T^4 
   *                                + a_9
   *    \f]
   *
   */

  template<typename CoeffType=double>
  class NASA9CurveFit
  {
  public:
    
    NASA9CurveFit( const std::vector<CoeffType>& coeffs );
    ~NASA9CurveFit();

    //! The number of intervals for this NASA9 curve fit
    unsigned int n_intervals() const;

    //! The interval the input temperature lies in
    /*!
      @returns which curve fit interval the input temperature 
      lies in.  The NASA9 thermodynamic intervals are 
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
            \frac{h}{\mathrm{R}T} = -\frac{a_0}{T^2} + \frac{a_1}{T} \ln(T) 
                                    + a2  + \frac{a_3}{2} T + \frac{a_4}{3} T^2 
                                    + \frac{a_5}{4} T^3 + \frac{a_6}{5} T^4 
                                    + \frac{a_8}{T}
        \f]
     */
    template <typename StateType>
    StateType h_over_RT( const TempCache<StateType>& cache) const;

    /*!
      @returns the value \f$\frac{s}{\mathrm{R}}\f$
        \f[
            \frac{s}{\mathrm{R}} = -\frac{a_0}{2T^2} - \frac{a_1}{T}
                                    + a_2 \ln(T) + a_3 T + \frac{a_4}{2} T^2 
                                    + \frac{a_5}{3} T^3 + \frac{a_6}{4} T^4 
                                    + a_9
        \f]
     */
    template <typename StateType>
    StateType s_over_R( const TempCache<StateType>& cache) const;

    /*!
      @returns the value \f$\frac{g}{\mathrm{R}T} = \frac{h}{\mathrm{R}T} - \frac{s}{R}\f$
        \f[
            \frac{g}{\mathrm{R}T} = -\frac{a_0}{2T^2} - \frac{a_1 + a_8}{T}
                                    + \frac{a_1}{T} \ln(T) - a2 \ln(T) 
                                    + (a_2 - a_9) - \frac{a_3}{2} T 
                                    - \frac{a_4}{6} T^2
                                    - \frac{a_5}{12} T^3 - \frac{a_6}{20} T^4 
        \f]
     */
    template <typename StateType>
    StateType h_RT_minus_s_R( const TempCache<StateType>& cache) const;

    /*!
      @returns the value \f$\frac{\partial\left(\frac{g}{\mathrm{R}T}\right)}\frac{\partial T} 
                            = \frac{\partial\left(\frac{h}{\mathrm{R}T} - \frac{s}{R}\right)}{\partial T}\f$
        \f[
            \frac{g}{\mathrm{R}T} =   \frac{a_0}{T^3} + \frac{a_1 + a_8}{T^2}
                                    + a_1 \ln(T) + a_1 - \frac{a2}{T}
                                    - \frac{a_3}{2} - \frac{a_4}{3} T
                                    - \frac{a_5}{4} T^2 - \frac{a_6}{5} T^3 
        \f]
     */
    template <typename StateType>
    StateType dh_RT_minus_s_R_dT( const TempCache<StateType>& cache) const;


    //! @returns a pointer to the coefficients in the interval specified.
    /*!   
      The NASA9-style equilibrium curve fits are defined in terms of
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

   private:
    // for compatibility with NASA, not used
    NASA9CurveFit( const std::vector<CoeffType>& coeffs, const std::vector<CoeffType> &temps );
  };


  /* ------------------------- Inline Functions -------------------------*/

  template<typename CoeffType>
  inline
  NASA9CurveFit<CoeffType>::NASA9CurveFit( const std::vector<CoeffType>& coeffs )
    : _n_coeffs(10),
      _coefficients(coeffs)
  {
    return;
  }


  template<typename CoeffType>
  inline
  NASA9CurveFit<CoeffType>::~NASA9CurveFit()
  {
    return;
  }


  template<typename CoeffType>
  template<typename StateType>
  inline
  typename Antioch::rebind<StateType, unsigned int>::type
  NASA9CurveFit<CoeffType>::interval(const StateType& T) const
  {
    typedef typename 
      Antioch::rebind<StateType, unsigned int>::type UIntType;
    UIntType interval;
    Antioch::zero_clone(interval, T);

    typedef typename Antioch::value_type<StateType>::type ScalarType;

    /* NASA9 thermodynamic intervals are:
       [200-1,000], [1,000-6,000], [6,000-20,000] K */
    interval = Antioch::if_else
      (T > ScalarType(6000.),
       Antioch::constant_clone(interval,2),
       UIntType
         (Antioch::if_else
            (T > ScalarType(1000.),
             Antioch::constant_clone(interval,1),
             Antioch::constant_clone(interval,0))));

    return interval;
  }


  template<typename CoeffType>
  inline
  unsigned int NASA9CurveFit<CoeffType>::n_intervals() const 
  { return _coefficients.size() / _n_coeffs; }


  template<typename CoeffType>
  inline
  const CoeffType* NASA9CurveFit<CoeffType>::coefficients(const unsigned int interval) const
  {
    antioch_assert_less( interval, this->n_intervals() );
    antioch_assert_less_equal( _n_coeffs*(interval+1), _coefficients.size() );
    
    return &_coefficients[_n_coeffs*interval];
  }

  template<typename CoeffType>
  template <typename StateType>
  inline
  const StateType NASA9CurveFit<CoeffType>::cp_over_R(const TempCache<StateType>& cache) const
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
           StateType(a[0]/cache.T2 + a[1]/cache.T + a[2] + a[3]*cache.T +
                     a[4]*cache.T2 + a[5]*cache.T3 + a[6]*cache.T4),
           returnval);
      }

    return returnval;
      
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType NASA9CurveFit<CoeffType>::h_over_RT( const TempCache<StateType>& cache) const
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
    
         /* h/RT = -a0*T^-2   + a1*T^-1*lnT + a2     + a3*T/2 + a4*T^2/3 + a5*T^3/4 + a6*T^4/5 + a8/T */
        returnval = Antioch::if_else
        ( interval == i,
           StateType( -a[0]/cache.T2 + a[1]*cache.lnT/cache.T + a[2] + 
                       a[3]*cache.T/2.0L + a[4]*cache.T2/3.0L + a[5]*cache.T3/4.0L +
                       a[6]*cache.T4/5.0L + a[8]/cache.T),
           returnval);
       }
       return returnval;
   }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType NASA9CurveFit<CoeffType>::s_over_R( const TempCache<StateType>& cache) const
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
    
    /* s/R = -a0*T^-2/2 - a1*T^-1     + a2*lnT + a3*T   + a4*T^2/2 + a5*T^3/3 + a6*T^4/4 + a9 */
        returnval = Antioch::if_else
        ( interval == i,
           StateType( -a[0]/cache.T2/2.0 - a[1]/cache.T + a[2]*cache.lnT 
                      + a[3]*cache.T + a[4]*cache.T2/2.0 + a[5]*cache.T3/3.0 
                      + a[6]*cache.T4/4.0 + a[9]),
           returnval);
       }
       return returnval;
   }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType
  NASA9CurveFit<CoeffType>::h_RT_minus_s_R( const TempCache<StateType>& cache) const
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
    
    /* h/RT = -a[0]/T2    + a[1]*lnT/T + a[2]     + a[3]*T/2. + a[4]*T2/3. + a[5]*T3/4. + a[6]*T4/5. + a[8]/T,
       s/R  = -a[0]/T2/2. - a[1]/T     + a[2]*lnT + a[3]*T    + a[4]*T2/2. + a[5]*T3/3. + a[6]*T4/4. + a[9]   */
        returnval = Antioch::if_else
        ( interval == i,
           StateType(-a[0]/cache.T2/2.0 + (a[1] + a[8])/cache.T +
                     a[1]*cache.lnT/cache.T - a[2]*cache.lnT + 
                     (a[2] - a[9]) - a[3]*cache.T/2.0 -
                     a[4]*cache.T2/6.0 - a[5]*cache.T3/12.0 -
                     a[6]*cache.T4/20.0),
           returnval);
       }
       return returnval;
   }

   template <typename CoeffType>
   template <typename StateType>
   inline
   StateType NASA9CurveFit<CoeffType>::dh_RT_minus_s_R_dT( const TempCache<StateType>& cache) const
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
           StateType(a[0]/cache.T3 - a[8]/cache.T2 -
                     a[1]*cache.lnT/cache.T2 - a[2]/cache.T -
                     a[3]/2.  - a[4]*cache.T/3. - a[5]*cache.T2/4. -
                     a[6]*cache.T3/5.),
           returnval);
      }

    return returnval;
        
   }


/*
    CEACurveFit is a synomym of NASA9CurveFit.

    Sylvain (06/02/2015 - european date format - dd/mm/yyyy):
    Food for the mind:
        the NASA7CurveFit has custom temperature
        ranges, while, as far as I know, the 
        9 coefficients form is solely encountered
        in the CEA program. This implementation
        acknowledge this and NASA9 and CEA are
        complete synomyms, there is no custom
        temperature ranges possible for them.

        As far as future developement goes, the
        question is ``will be want someday to
        pick that lock?''. The thing at stake here is
        the possible antioch_deprecated() in the
        CEA constructor. As of now, they're full
        synomyms, so an antioch_deprecated() should
        be in order, but I personnaly prefers more
        general things, specially if we can sort
        them out at compile time by template and
        metaprogramming shims. I'm all for it Paul
        is not, history shows I win if I find a
        real case where it's needed...

        So I'm choosing the middle ground: no
        antioch_deprecated() when they are
        full synomyms.
*/


  /*!\class CEACurveFit

     This is a synomym for NASA9CurveFit, the NASA9
     name ensures consistency with the NASA7 objects
     while the CEA name provides backward compatiblity
     and a more `physics-based' name.

     Note that as nothing happens in the destructor (and
     nothing should ever), no need to get virtual in NASA9.
  */
  template<typename CoeffType=double>
  class CEACurveFit:public NASA9CurveFit<CoeffType>
  {
      CEACurveFit( const std::vector<CoeffType>& coeffs ):NASA9CurveFit(coeffs){}
      ~CEACurveFit(){}
  }

} // end namespace Antioch

#endif //ANTIOCH_NASA9_CURVE_FIT_H
