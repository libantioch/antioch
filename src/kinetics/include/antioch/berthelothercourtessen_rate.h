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

// C++
#include <cmath>
#include <iostream>
#include <sstream>

namespace Antioch
{
  //! Berthelot Hercourt-Essen rate equation.
  /*!
   * Berthelot Hercourt-Essen rate equation.  Computes rates of the form
   * \f$ C_f\times T^\eta\times \exp(D*T) \f$.
   */
  template<typename CoeffType=double>
  class BerthelotHercourtEssenRate: public KineticsType<CoeffType>
  {
  
  public:

    BerthelotHercourtEssenRate (const CoeffType Cf=0., const CoeffType eta=0., const CoeffType D=0., const CoeffType Tref = 1.);
    ~BerthelotHercourtEssenRate();
    
    void set_Cf(  const CoeffType Cf );
    void set_eta( const CoeffType eta );
    void set_D(   const CoeffType D );
    void set_Tref(const CoeffType Tref );

    CoeffType Cf()   const;
    CoeffType eta()  const;
    CoeffType D()    const;
    CoeffType Tref() const;

    //! \return the rate evaluated at \p T.
    template <typename StateType>
    StateType operator()(const StateType& T) const;

    //! \return the rate evaluated at \p T.
    template <typename StateType>
    StateType rate(const StateType& T) const;

    //! \return the derivative with respect to temperature evaluated at \p T.
    template <typename StateType>
    StateType derivative( const StateType& T ) const;

    //! Simultaneously evaluate the rate and its derivative at \p T.
    template <typename StateType>
    void rate_and_derivative(const StateType& T, StateType& rate, StateType& drate_dT) const;

    //! print equation
    const std::string numeric() const;

  private:

    void compute_cf();

    CoeffType _raw_Cf;
    CoeffType _Cf;
    CoeffType _eta;
    CoeffType _D;
    CoeffType _Tref;
    
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
  StateType BerthelotHercourtEssenRate<CoeffType>::operator()(const StateType& T) const
  {
    using std::pow;
    using std::exp;
    return _Cf* (pow(T,_eta)*exp(_D*T));
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType BerthelotHercourtEssenRate<CoeffType>::rate(const StateType& T) const
  {
    using std::pow;
    using std::exp;
    return _Cf* (pow(T,_eta)*exp(_D*T));
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType BerthelotHercourtEssenRate<CoeffType>::derivative( const StateType& T ) const
  {
    return (*this)(T)*(_eta/T + _D);
  }

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
