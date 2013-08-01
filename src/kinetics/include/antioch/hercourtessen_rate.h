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

#ifndef ANTIOCH_HERCOURT_ESSEN_RATE_H
#define ANTIOCH_HERCOURT_ESSEN_RATE_H

//Antioch
#include "antioch/kinetics_type.h"
#include "antioch/cmath_shims.h"

// C++
#include <cmath>
#include <iostream>
#include <sstream>

namespace Antioch
{
  //! Hercourt-Essen rate equation.
  /*!
   * Hercourt-Essen rate equation.  Computes rates of the form
   * \f$ C_f\times T^\eta \f$.
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
    _Cf = _raw_Cf * pow(KineticsModel::Tref<CoeffType>()/_Tref,_eta);

    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_HERCOURT_ESSEN_RATE_H
