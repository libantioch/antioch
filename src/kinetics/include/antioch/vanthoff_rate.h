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

#ifndef ANTIOCH_VAN_T_HOFF_RATE_H
#define ANTIOCH_VAN_T_HOFF_RATE_H

//Antioch
#include "antioch/kinetics_type.h"
#include "antioch/physical_constants.h"

// C++
#include <cmath>
#include <iostream>
#include <sstream>

namespace Antioch
{
  //! Van't Hoff rate equation.
  /*!
   * 
   * Van't Hoff rate equation.  Computes rates of the form
   * \f$ C_f\times T^\eta\times \exp(-E_a/T + DT) \f$.
   */
  template<typename CoeffType=double>
  class VantHoffRate:public KineticsType<CoeffType>
  {
  
  public:

    VantHoffRate (const CoeffType Cf=0., const CoeffType eta=0., const CoeffType Ea=0., const CoeffType D=0., const CoeffType Tref = 1., const CoeffType rscale = Constants::R_universal<CoeffType>()/1000.);
    ~VantHoffRate();
    
    void set_Cf(    const CoeffType Cf );
    void set_eta(   const CoeffType eta );
    void set_Ea(    const CoeffType Ea );
    void set_D(     const CoeffType D );
    void set_Tref(  const CoeffType Tref );
    void set_rscale(const CoeffType rscale );

    CoeffType Cf()     const;
    CoeffType eta()    const;
    CoeffType Ea()     const;
    CoeffType D()      const;
    CoeffType Tref()   const;
    CoeffType rscale() const;

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

    CoeffType _raw_Cf;
    CoeffType _Cf;
    CoeffType _eta;
    CoeffType _raw_Ea;
    CoeffType _Ea;
    CoeffType _D;
    CoeffType _Tref;
    CoeffType _rscale;
    
  };

  template<typename CoeffType>
  VantHoffRate<CoeffType>::VantHoffRate(const CoeffType Cf, const CoeffType eta, const CoeffType Ea, const CoeffType D, const CoeffType Tref, const CoeffType rscale)
    : KineticsType<CoeffType>(KineticsModel::VANTHOFF),
      _raw_Cf(Cf),
      _eta(eta),
      _raw_Ea(Ea),
      _D(D),
      _Tref(Tref),
      _rscale(rscale)
  {
    using std::pow;
    _Ea = _raw_Ea / _rscale;
    _Cf = _raw_Cf * pow(KineticsModel::Tref/_Tref,_eta);
    return;
  }

  template<typename CoeffType>
  VantHoffRate<CoeffType>::~VantHoffRate()
  {
    return;
  }

  template<typename CoeffType>
  const std::string VantHoffRate<CoeffType>::numeric() const
  {
    std::stringstream os;
    os << _raw_Cf;
    if (_eta != 0.) os << "*(T/" << _Tref << ")^" << _eta;
    os << "*exp(-" << _raw_Ea << "/(R*T) + " << _D << "*T)";

    return os.str();
  }

  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType>
  inline
  void VantHoffRate<CoeffType>::set_Cf( const CoeffType Cf )
  {
    using std::pow;
    _raw_Cf = Cf;
    _Cf = _raw_Cf * pow(KineticsModel::Tref/_Tref,_eta);
    return;
  }

  template<typename CoeffType>
  inline
  void VantHoffRate<CoeffType>::set_Tref( const CoeffType Tref )
  {
    using std::pow;
    _Tref = Tref;
    _Cf = _raw_Cf * pow(KineticsModel::Tref/_Tref,_eta);
    return;
  }

  template<typename CoeffType>
  inline
  void VantHoffRate<CoeffType>::set_eta( const CoeffType eta )
  {
    _eta = eta;
    return;
  }

  template<typename CoeffType>
  inline
  void VantHoffRate<CoeffType>::set_Ea( const CoeffType Ea )
  {
    _raw_Ea = Ea;
    _Ea = _raw_Ea / _rscale;
    return;
  }

  template<typename CoeffType>
  inline
  void VantHoffRate<CoeffType>::set_rscale( const CoeffType rscale )
  {
    _rscale = rscale;
    _Ea = _raw_Ea / _rscale;
    return;
  }

  template<typename CoeffType>
  inline
  void VantHoffRate<CoeffType>::set_D( const CoeffType D )
  {
    _D = D;
    return;
  }


  template<typename CoeffType>
  inline
  CoeffType VantHoffRate<CoeffType>::Cf() const
  { return _Cf; }

  template<typename CoeffType>
  inline
  CoeffType VantHoffRate<CoeffType>::eta() const
  { return _eta; }

  template<typename CoeffType>
  inline
  CoeffType VantHoffRate<CoeffType>::Ea() const
  { return _Ea; }

  template<typename CoeffType>
  inline
  CoeffType VantHoffRate<CoeffType>::D() const
  { return _D; }

  template<typename CoeffType>
  inline
  CoeffType VantHoffRate<CoeffType>::rscale() const
  { return _rscale; }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType VantHoffRate<CoeffType>::operator()(const StateType& T) const
  {
    using std::pow;
    using std::exp;
    return _Cf* (pow(T,_eta)*exp(-_Ea/T + _D*T));
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType VantHoffRate<CoeffType>::rate(const StateType& T) const
  {
    using std::pow;
    using std::exp;
    return _Cf* (pow(T,_eta)*exp(-_Ea/T + _D*T));
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType VantHoffRate<CoeffType>::derivative( const StateType& T ) const
  {
    return (*this)(T)*(_D + _eta/T + _Ea/(T*T));
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  void VantHoffRate<CoeffType>::rate_and_derivative( const StateType& T,
						      StateType& rate,
						      StateType& drate_dT) const
  {
    rate     = (*this)(T);
    drate_dT = rate*(_D + _eta/T + _Ea/(T*T));
    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_VAN_T_HOFF_RATE_H
