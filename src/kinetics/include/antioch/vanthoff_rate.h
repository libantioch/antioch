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

// C++
#include <cmath>
#include <iostream>

namespace Antioch
{
  //! Van't Hoff rate equation.
  /*!
   * 
   * Van't Hoff rate equation.  Computes rates of the form
   * \f$ C_f\times T^\eta\times \exp(-E_a/T + DT) \f$.
   */
  template<typename CoeffType=double>
  class VantHoffRate:public KineticsType
  {
  
  public:

    VantHoffRate (const CoeffType Cf=0., const CoeffType eta=0., const CoeffType Ea=0., const CoeffType D=0.);
    ~VantHoffRate();
    
    void set_Cf( const CoeffType Cf );
    void set_eta( const CoeffType eta );
    void set_Ea( const CoeffType Ea );
    void set_D( const CoeffType D );

    void scale_Ea( const CoeffType scale );
    void scale_D( const CoeffType scale );

    CoeffType Cf()  const;
    CoeffType eta() const;
    CoeffType Ea()  const;
    CoeffType D()   const;

    //! \return the rate evaluated at \p T.
    template <typename StateType>
    StateType operator()(const StateType& T) const;

    //! \return the derivative with respect to temperature evaluated at \p T.
    template <typename StateType>
    StateType derivative( const StateType& T ) const;

    //! Simultaneously evaluate the rate and its derivative at \p T.
    template <typename StateType>
    void rate_and_derivative(const StateType& T, StateType& rate, StateType& drate_dT) const;

    //! Formatted print, by default to \p std::cout
    void print(std::ostream& os = std::cout) const;

    //! Formatted print.
    friend std::ostream& operator<<(std::ostream& os, const VantHoffRate& rate)
    {
      rate.print(os);
      return os;
    }

  private:

    CoeffType _Cf;
    CoeffType _eta;
    CoeffType _Ea;
    CoeffType _D;
    
  };

  template<typename CoeffType>
  VantHoffRate<CoeffType>::VantHoffRate(const CoeffType Cf, const CoeffType eta, const CoeffType Ea, const CoeffType D)
    : _Cf(Cf),
      _eta(eta),
      _Ea(Ea),
      _D(D)
  {
    return;
  }

  template<typename CoeffType>
  VantHoffRate<CoeffType>::~VantHoffRate()
  {
    return;
  }

  template<typename CoeffType>
  void VantHoffRate<CoeffType>::print(std::ostream& os) const
  {
    os << _Cf;
    if (_eta != 0.) os << "*T^" << _eta;
    os << "*exp(-" << _Ea << "/T + " << _D << "*T)";

    return;
  }

  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType>
  inline
  void VantHoffRate<CoeffType>::set_Cf( const CoeffType Cf )
  {
    _Cf = Cf;
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
    _Ea = Ea;
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
  void VantHoffRate<CoeffType>::scale_Ea( const CoeffType scale )
  {
    _Ea *= scale;
    return;
  }

  template<typename CoeffType>
  inline
  void VantHoffRate<CoeffType>::scale_D( const CoeffType scale )
  {
    _D *= scale;
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

#endif // ANTIOCH_VAN-T-HOFF_RATE_H
