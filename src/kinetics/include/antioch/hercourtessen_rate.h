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

#ifndef ANTIOCH_HERCOURT_ESSEN_RATE_H
#define ANTIOCH_HERCOURT_ESSEN_RATE_H

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/kinetics_type.h"
#include "antioch/cmath_shims.h"

// C++
#include <cmath>
#include <iostream>
#include <sstream>

namespace Antioch
{
  //! Hercourt-Essen rate equation.
  /*!\class HercourtEssenRate
 *
   * Hercourt-Essen is a kinetics model of the form:
   * \f[
   *  \alpha(T) = A \left(\frac{T}{\mathrm{T_\text{ref}}}\right)^\beta
   * \f]
   * Thus
   * \f[
   *  \frac{\partial \alpha(T)}{\partial T} = \alpha(T) \frac{\beta}{T}
   * \f]
   *
   * Internally, we use the reduced pre-exponential factor \f$a = \frac{A}{\mathrm{T_\text{ref}}^{\beta}}\f$
   */
  template<typename CoeffType=double>
  class HercourtEssenRate : public KineticsType<CoeffType>
  {

  private:
  
    CoeffType _raw_Cf;
    CoeffType _Cf;
    CoeffType _eta;
    CoeffType _Tref;
    
  public:

    HercourtEssenRate (const CoeffType Cf=0., const CoeffType eta=0., const CoeffType Tref = 1.);
    ~HercourtEssenRate();
    
    void set_Cf(  const CoeffType Cf );
    void set_eta( const CoeffType eta );
    void set_Tref(const CoeffType Tref );

    /*! reset the coeffs
     *
     * Two ways of modifying your rate:
     *   - you change totally the rate, thus you 
     *        require exactly three parameters, the order
     *        assumed is Cf, eta, Tref
     *   - you just change the value, thus Tref is not
     *        modified. You require exactly two parameters,
     *        the order assumed is Cf, eta
     */
    template <typename VectorCoeffType>
    void reset_coefs(const VectorCoeffType & coefficients);

    CoeffType Cf()   const;
    CoeffType eta()  const;
    CoeffType Tref() const;

    //! \return the rate evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    rate(const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, _Cf* (ant_pow(T,_eta)))

    //! \return the rate evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    operator()(const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, this->rate(T))

    //! \return the derivative with respect to temperature evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    derivative( const StateType& T ) const
    ANTIOCH_AUTOFUNC(StateType, (*this)(T)/T*(_eta))

    //! Simultaneously evaluate the rate and its derivative at \p T.
    template <typename StateType>
    void rate_and_derivative(const StateType& T, StateType& rate, StateType& drate_dT) const;

    //! print equation
    const std::string numeric() const;

  private:

    void compute_cf();

  };

  template<typename CoeffType>
  HercourtEssenRate<CoeffType>::HercourtEssenRate(const CoeffType Cf, const CoeffType eta, const CoeffType Tref)
    : KineticsType<CoeffType>(KineticsModel::HERCOURT_ESSEN),
      _raw_Cf(Cf),
      _eta(eta),
      _Tref(Tref)
  {
    this->compute_cf();

    return;
  }

  template<typename CoeffType>
  HercourtEssenRate<CoeffType>::~HercourtEssenRate()
  {
    return;
  }

  template<typename CoeffType>
  const std::string HercourtEssenRate<CoeffType>::numeric() const
  {
    std::stringstream os;
    os << _raw_Cf;
    os << "*(T/" << _Tref << ")^" << _eta;

    return os.str();
  }

  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType>
  inline
  void HercourtEssenRate<CoeffType>::set_Cf( const CoeffType Cf )
  {
    _raw_Cf = Cf;
    this->compute_cf();

    return;
  }

  template<typename CoeffType>
  inline
  void HercourtEssenRate<CoeffType>::set_eta( const CoeffType eta )
  {
    _eta = eta;
    this->compute_cf();
    return;
  }

  template<typename CoeffType>
  inline
  void HercourtEssenRate<CoeffType>::set_Tref( const CoeffType Tref )
  {
    _Tref = Tref;
    this->compute_cf();

    return;
  }

  template<typename CoeffType>
  template <typename VectorCoeffType>
  inline
  void HercourtEssenRate<CoeffType>::reset_coefs(const VectorCoeffType & coefficients)
  {
      antioch_assert_greater(coefficients.size(),1);
      antioch_assert_less(coefficients.size(),4);
      if(coefficients.size() == 3)this->set_Tref(coefficients[2]);
      this->set_Cf(coefficients[0]);
      this->set_eta(coefficients[1]);
  }

  template<typename CoeffType>
  inline
  CoeffType HercourtEssenRate<CoeffType>::Cf() const
  { return _Cf; }

  template<typename CoeffType>
  inline
  CoeffType HercourtEssenRate<CoeffType>::eta() const
  { return _eta; }

  template<typename CoeffType>
  inline
  CoeffType HercourtEssenRate<CoeffType>::Tref() const
  { return _Tref; }

  template<typename CoeffType>
  template<typename StateType>
  inline
  void HercourtEssenRate<CoeffType>::rate_and_derivative( const StateType& T,
                                                          StateType& rate,
                                                          StateType& drate_dT) const
  {
    rate     = (*this)(T);
    drate_dT = rate/T*(_eta);
    return;
  }

  template<typename CoeffType>
  inline
  void HercourtEssenRate<CoeffType>::compute_cf()
  {
    _Cf = _raw_Cf * ant_pow(KineticsModel::Tref<CoeffType>()/_Tref,_eta);

    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_HERCOURT_ESSEN_RATE_H
