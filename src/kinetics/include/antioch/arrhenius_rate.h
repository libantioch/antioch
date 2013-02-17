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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef ANTIOCH_ARRHENIUS_RATE_H
#define ANTIOCH_ARRHENIUS_RATE_H

// C++
#include <cmath>
#include <iostream>

namespace Antioch
{
  //! Modified Arrhenius rate equation.
  /*!
   * Modified Arrhenius rate equation.  Computes rates of the form
   * \f$ C_f\times T^\eta\times \exp(-E_a/T) \f$. This class copied from
   * the \p FIN-S code and slightly reformatted for \p Antioch.
   */
  template<class NumericType>
  class ArrheniusRate
  {
  
  public:

    ArrheniusRate (const NumericType Cf=0., const NumericType eta=0., const NumericType Ea=0.);
    ~ArrheniusRate();
    
    void set_Cf( const NumericType Cf );
    void set_eta( const NumericType eta );
    void set_Ea( const NumericType Ea );

    void scale_Ea( const NumericType scale );

    NumericType Cf() const;
    NumericType eta() const;
    NumericType Ea() const;

    //! \return the rate evaluated at \p T.
    NumericType operator()(const NumericType T) const;

    //! \return the derivative with respect to temperature evaluated at \p T.
    NumericType derivative( const NumericType T ) const;

    //! Simultaneously evaluate the rate and its derivative at \p T.
    void rate_and_derivative(const NumericType T, NumericType& rate, NumericType& drate_dT) const;

    //! Formatted print, by default to \p std::cout
    void print(std::ostream& os = std::cout) const;

    //! Formatted print.
    friend std::ostream& operator<<(std::ostream& os, const ArrheniusRate& rate)
    {
      rate.print(os);
      return os;
    }

  private:

    NumericType _Cf;
    NumericType _eta;
    NumericType _Ea;
    
  };

  template<class NumericType>
  ArrheniusRate<NumericType>::ArrheniusRate(const NumericType Cf, const NumericType eta, const NumericType Ea)
    : _Cf(Cf),
      _eta(eta),
      _Ea(Ea)
  {
    return;
  }

  template<class NumericType>
  ArrheniusRate<NumericType>::~ArrheniusRate()
  {
    return;
  }

  template<class NumericType>
  void ArrheniusRate<NumericType>::print(std::ostream& os) const
  {
    os << _Cf;
    if (_eta != 0.) os << "*T^" << _eta;
    os << "*exp(-" << _Ea << "/T)";

    return;
  }

  /* ------------------------- Inline Functions -------------------------*/
  template<class NumericType>
  inline
  void ArrheniusRate<NumericType>::set_Cf( const NumericType Cf )
  {
    _Cf = Cf;
    return;
  }

  template<class NumericType>
  inline
  void ArrheniusRate<NumericType>::set_eta( const NumericType eta )
  {
    _eta = eta;
    return;
  }

  template<class NumericType>
  inline
  void ArrheniusRate<NumericType>::set_Ea( const NumericType Ea )
  {
    _Ea = Ea;
    return;
  }

  template<class NumericType>
  inline
  void ArrheniusRate<NumericType>::scale_Ea( const NumericType scale )
  {
    _Ea *= scale;
    return;
  }

  template<class NumericType>
  inline
  NumericType ArrheniusRate<NumericType>::Cf() const
  { return _Cf; }

  template<class NumericType>
  inline
  NumericType ArrheniusRate<NumericType>::eta() const
  { return _eta; }

  template<class NumericType>
  inline
  NumericType ArrheniusRate<NumericType>::Ea() const
  { return _Ea; }

  template<class NumericType>
  inline
  NumericType ArrheniusRate<NumericType>::operator()(const NumericType T) const
  {
    return _Cf* (std::pow(T,_eta)*std::exp(-_Ea/T));
  }

  template<class NumericType>
  inline
  NumericType ArrheniusRate<NumericType>::derivative( const NumericType T ) const
  {
    return (*this)(T)/T*(_eta + _Ea/T);
  }

  template<class NumericType>
  inline
  void ArrheniusRate<NumericType>::rate_and_derivative( const NumericType T,
							NumericType& rate,
							NumericType& drate_dT) const
  {
    rate     = (*this)(T);
    drate_dT = rate/T*(_eta + _Ea/T);
    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_ARRHENIUS_RATE_H
