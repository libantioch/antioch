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

#ifndef ANTIOCH_BERTHELOT_HERCOURT_ESSEN_RATE_H
#define ANTIOCH_BERTHELOT_HERCOURT_ESSEN_RATE_H

//Antioch
#include "antioch/kinetics_type.h"
#include "antioch/cmath_shims.h"

// C++
#include <cmath>
#include <iostream>
#include <sstream>

namespace Antioch
{
  //! Berthelot Hercourt-Essen rate equation.
  /*!\class BerthelotHercourtEssenRate
 *
   * The Berthelot Hercourt Essen kinetics model is of the form:
   * \f[
   *    \alpha(T) =  C_f \left(\frac{T}{\mathrm{T_\text{ref}}}\right)^\beta  \exp\left(D T\right) 
   * \f]
   * thus
   * \f[
   *    \frac{\partial \alpha(T)}{\partial T} = \alpha(T) \left(\frac{\beta}{T} + D\right)
   * \f]
   *
   * Internally, we use the reduced temperature \f$t = \frac{T}{\mathrm{T_\text{ref}}}\f$.
   */
  template<typename CoeffType=double>
  class BerthelotHercourtEssenRate: public KineticsType<CoeffType>
  {

  private:
  
    CoeffType _raw_Cf;
    CoeffType _Cf;
    CoeffType _eta;
    CoeffType _D;
    CoeffType _Tref;
 
  public:

    BerthelotHercourtEssenRate (const CoeffType Cf=0., const CoeffType eta=0., const CoeffType D=0., const CoeffType Tref = 1.);
    ~BerthelotHercourtEssenRate();
   
    /** \brief set the pre-exponential factor in units {kmol,s,m} (e.g. order=1 : [s^-1]; order=2 : [m^3.kmol^-1.s^-1] )
     * 
     * \param Cf  pre-exponential factor in units {kmol,s,m} (e.g. order=1 : [s^-1]; order=2 : [m^3.kmol^-1.s^-1] )
     */ 
    void set_Cf(  const CoeffType Cf );
    
     /** \brief set \f$ \beta \f$ of \f$ \left(\frac{T}{\mathrm{T_\text{ref}}}\right)^\beta \f$
     * 
     * \param eta \f$ \beta \f$ value [none] 
     */  
    void set_eta( const CoeffType eta );
    
    /** \brief set the exponential factor 
     * 
     * \param D  factor in [K^-1]
     */ 
    void set_D(   const CoeffType D );
    
     /** \brief set the reference temperature
     * 
     * \param Tref  reference temperature in [K]
     */    
    void set_Tref(const CoeffType Tref );

    //! \return the pre-exponential factor in units {kmol,s,m} (e.g. order=1 : [s^-1]; order=2 : [m^3.kmol^-1.s^-1] )
    CoeffType Cf()   const;
    
    //! \return the  \f$ \beta \f$ of \f$ \left(\frac{T}{\mathrm{T_\text{ref}}}\right)^\beta \f$ in [none]
    CoeffType eta()  const;
    
    //! \return the exponential factor in [K^-1]
    CoeffType D()    const;
    
    //! \return reference temperature in [K]
    CoeffType Tref() const;

    //! \return the rate evaluated at \p T in units {kmol,s,m} (e.g. order=1 : [s^-1]; order=2 : [m^3.kmol^-1.s^-1] ).
    template <typename StateType>
    ANTIOCH_AUTO(StateType)
    rate(const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, _Cf*(ant_pow(T,_eta)*ant_exp(_D*T)))

    //! \return the rate evaluated at \p T in units {kmol,s,m} (e.g. order=1 : [s^-1]; order=2 : [m^3.kmol^-1.s^-1] ).
    template <typename StateType>
    ANTIOCH_AUTO(StateType)
    operator()(const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, this->rate(T))

    //! \return the derivative with respect to temperature evaluated at \p T in units {kmol,s,m,K}.
    template <typename StateType>
    ANTIOCH_AUTO(StateType)
    derivative( const StateType& T ) const
    ANTIOCH_AUTOFUNC(StateType, (*this)(T)*(_eta/T + _D))

    //! Simultaneously evaluate the rate and its derivative at \p T in units {kmol,s,m,K}.
    template <typename StateType>
    void rate_and_derivative(const StateType& T, StateType& rate, StateType& drate_dT) const;

    //! print equation
    const std::string numeric() const;

  private:

    void compute_cf();
  };

  template<typename CoeffType>
  BerthelotHercourtEssenRate<CoeffType>::BerthelotHercourtEssenRate(const CoeffType Cf, const CoeffType eta, const CoeffType D, const CoeffType Tref)
    : KineticsType<CoeffType>(KineticsModel::BHE),
      _raw_Cf(Cf),
      _eta(eta),
      _D(D),
      _Tref(Tref)
  {
    using std::pow;
    this->compute_cf();
    return;
  }

  template<typename CoeffType>
  BerthelotHercourtEssenRate<CoeffType>::~BerthelotHercourtEssenRate()
  {
    return;
  }

  template<typename CoeffType>
  const std::string BerthelotHercourtEssenRate<CoeffType>::numeric() const
  {
    std::stringstream os;
    os << _raw_Cf;
    if (_eta != 0.) os << "*(T/" << _Tref << ")^" << _eta;
    os << "*exp(" << _D << "*T)";

    return os.str();
  }

  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType>
  inline
  void BerthelotHercourtEssenRate<CoeffType>::set_Cf( const CoeffType Cf )
  {
    using std::pow;
    _raw_Cf = Cf;
    this->compute_cf();
    return;
  }

  template<typename CoeffType>
  inline
  void BerthelotHercourtEssenRate<CoeffType>::set_Tref( const CoeffType Tref )
  {
    using std::pow;
    _Tref = Tref;
    this->compute_cf();
    return;
  }

  template<typename CoeffType>
  inline
  void BerthelotHercourtEssenRate<CoeffType>::set_eta( const CoeffType eta )
  {
    _eta = eta;
    return;
  }

  template<typename CoeffType>
  inline
  void BerthelotHercourtEssenRate<CoeffType>::set_D( const CoeffType D )
  {
    _D = D;
    return;
  }

  template<typename CoeffType>
  inline
  CoeffType BerthelotHercourtEssenRate<CoeffType>::Cf() const
  { return _Cf; }

  template<typename CoeffType>
  inline
  CoeffType BerthelotHercourtEssenRate<CoeffType>::eta() const
  { return _eta; }

  template<typename CoeffType>
  inline
  CoeffType BerthelotHercourtEssenRate<CoeffType>::D() const
  { return _D; }

  template<typename CoeffType>
  inline
  CoeffType BerthelotHercourtEssenRate<CoeffType>::Tref() const
  { return _Tref; }

  template<typename CoeffType>
  template<typename StateType>
  inline
  void BerthelotHercourtEssenRate<CoeffType>::rate_and_derivative( const StateType& T,
                                                                   StateType& rate,
                                                                   StateType& drate_dT) const
  {
    rate     = (*this)(T);
    drate_dT = rate*(_eta/T + _D);
    return;
  }

  template<typename CoeffType>
  inline
  void BerthelotHercourtEssenRate<CoeffType>::compute_cf()
  {
    _Cf = _raw_Cf * pow(KineticsModel::Tref<CoeffType>()/_Tref,_eta);
    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_BERTHELOT_HERCOURT_ESSEN_RATE_H
