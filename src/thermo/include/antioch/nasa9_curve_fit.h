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

#ifndef ANTIOCH_NASA9_CURVE_FIT_H
#define ANTIOCH_NASA9_CURVE_FIT_H

// Antioch
#include "antioch/nasa_curve_fit_base.h"
#include "antioch/temp_cache.h"

namespace Antioch
{
  /*! \class NASA9CurveFit
   *
   *  This class stores the NASA polynomial fit to
   *  the thermodynamics quantities \f$\frac{C_p}{\mathrm{R}}\f$
   *  \f$\frac{h}{\mathrm{R}T}\f$ and \f$\frac{s}{\mathrm{R}T}\f$.
   *  This formulation requires nine coefficients, from \f$a_0\f$
   *  to \f$a_8\f$.
   *
   *  The default temperature intervals are
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
   *                                + \frac{a_7}{T}
   *    \f]
   *
   *    \f[
   *        \frac{s}{\mathrm{R}} = -\frac{a_0}{2T^2} - \frac{a_1}{T}
   *                                + a_2 \ln(T) + a_3 T + \frac{a_4}{2} T^2
   *                                + \frac{a_5}{3} T^3 + \frac{a_6}{4} T^4
   *                                + a_8
   *    \f]
   *
   */

  template<typename CoeffType=double>
  class NASA9CurveFit : public NASACurveFitBase<CoeffType>
  {
  public:

    //! Accepts a vector temperatures to specify the temp. intervals
    NASA9CurveFit( const std::vector<CoeffType>& coeffs, const std::vector<CoeffType>& temps );

    //! Assumes the coefficients correspond to the default temperature intervals
    /*! If there are not enough coefficients for the third interval,
        then this only builds up through the second interval. */
    NASA9CurveFit( const std::vector<CoeffType>& coeffs );

    ~NASA9CurveFit(){};

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
      @returns the derivative of \f$\frac{Cp}{\mathrm{R}}\f$ with
      respect to temperature
        \f[
           \frac{Cp}{\mathrm{R}} = \frac{-2 a_0}{T^3} - \frac{a_1}{T^2} + a_3 +
                                    2 a_4 T + 3 a_5 T^2 + 4 a_6 T^3
        \f]
     */
    template <typename StateType>
    const StateType dcp_over_R_dT (const TempCache<StateType> & cache) const;

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
      @returns the value \f$\frac{\partial\left(\frac{g}{\mathrm{R}T}\right)}{\partial T}
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

  protected:

    NASA9CurveFit(){};

    void init_nasa9_temps( const std::vector<CoeffType>& coeffs,
                           unsigned n_coeffs );

  };


  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType>
  inline
  NASA9CurveFit<CoeffType>::NASA9CurveFit( const std::vector<CoeffType>& coeffs,
                                           const std::vector<CoeffType>& temp )
    : NASACurveFitBase<CoeffType>(coeffs,temp)
  {
    this->_n_coeffs = 9;

    this->check_coeff_size();
    this->check_temp_coeff_size_consistency();
  }

  template<typename CoeffType>
  inline
  NASA9CurveFit<CoeffType>::NASA9CurveFit( const std::vector<CoeffType>& coeffs )
    : NASACurveFitBase<CoeffType>(coeffs,std::vector<CoeffType>())
  {
    this->_n_coeffs = 9;
    this->check_coeff_size();

    this->init_nasa9_temps( coeffs, this->_n_coeffs );

    this->check_temp_coeff_size_consistency();
  }

  template<typename CoeffType>
  inline
  void NASA9CurveFit<CoeffType>::init_nasa9_temps( const std::vector<CoeffType>& coeffs,
                                                   unsigned n_coeffs )
  {
    // All default NASA9 curve fits should have this
    this->_temp.resize(3);
    this->_temp[0] = 200.L;
    this->_temp[1] = 1000.L;
    this->_temp[2] = 6000.L;

    // But some also have this last interval. We infer this from the size of the coeffs.
    if( coeffs.size()/n_coeffs == 3 )
      this->_temp.push_back(20000.L);
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
  template <typename StateType>
  inline
  const StateType NASA9CurveFit<CoeffType>::dcp_over_R_dT(const TempCache<StateType>& cache) const
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
           StateType(-2*a[0]/cache.T3 - a[1]/cache.T2 + a[3] +
                     2*a[4]*cache.T + 3*a[5]*cache.T2 + 4*a[6]*cache.T3),
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
         const CoeffType *a = this->coefficients(i);

         /* h/RT = -a0*T^-2   + a1*T^-1*lnT + a2     + a3*T/2 + a4*T^2/3 + a5*T^3/4 + a6*T^4/5 + a7/T */
        returnval = Antioch::if_else
        ( interval == i,
           StateType( -a[0]/cache.T2 +
                      a[1]*cache.lnT/cache.T +
                      a[2] +
                      a[3]*cache.T/2 +
                      a[4]*cache.T2/3 +
                      a[5]*cache.T3/4 +
                      a[6]*cache.T4/5 +
                      a[7]/cache.T ),
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
         const CoeffType *a = this->coefficients(i);

    /* s/R = -a0*T^-2/2 - a1*T^-1     + a2*lnT + a3*T   + a4*T^2/2 + a5*T^3/3 + a6*T^4/4 + a8 */
        returnval = Antioch::if_else
        ( interval == i,
          StateType( -a[0]/cache.T2/2 -
                     a[1]/cache.T +
                     a[2]*cache.lnT +
                     a[3]*cache.T +
                     a[4]*cache.T2/2 +
                     a[5]*cache.T3/3 +
                     a[6]*cache.T4/4 +
                     a[8] ),
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
         const CoeffType *a = this->coefficients(i);

    /* h/RT = -a[0]/T2    + a[1]*lnT/T + a[2]     + a[3]*T/2. + a[4]*T2/3. + a[5]*T3/4. + a[6]*T4/5. + a[7]/T,
       s/R  = -a[0]/T2/2. - a[1]/T     + a[2]*lnT + a[3]*T    + a[4]*T2/2. + a[5]*T3/3. + a[6]*T4/4. + a[8]   */
        returnval = Antioch::if_else
        ( interval == i,
          StateType(-a[0]/cache.T2/2
                    + (a[1] + a[7])/cache.T
                    + a[1]*cache.lnT/cache.T
                    - a[2]*cache.lnT
                    + (a[2] - a[8])
                    - a[3]*cache.T/2
                    - a[4]*cache.T2/6
                    - a[5]*cache.T3/12
                    - a[6]*cache.T4/20 ),
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

    /* h/RT = -a[0]/T2    + a[1]*lnT/T + a[2]     + a[3]*T/2. + a[4]*T2/3. + a[5]*T3/4. + a[6]*T4/5. + a[7]/T,
       s/R  = -a[0]/T2/2. - a[1]/T     + a[2]*lnT + a[3]*T    + a[4]*T2/2. + a[5]*T3/3. + a[6]*T4/4. + a[8]   */
    for (unsigned int i=begin_interval; i != end_interval; ++i)
      {
        const CoeffType * const a =
          this->coefficients(i);
        returnval = Antioch::if_else
          (interval == i,
           StateType(a[0]/cache.T3 - a[7]/cache.T2 -
                     a[1]*cache.lnT/cache.T2 - a[2]/cache.T -
                     a[3]/2  - a[4]*cache.T/3 - a[5]*cache.T2/4 -
                     a[6]*cache.T3/5 ),
           returnval);
      }

    return returnval;

   }

} // end namespace Antioch
#endif // ANTIOCH_NASA9_CURVE_FIT_H
