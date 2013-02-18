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
  template<typename CoeffType=double>
  class ArrheniusRate
  {
  
  public:

    ArrheniusRate (const CoeffType Cf=0., const CoeffType eta=0., const CoeffType Ea=0.);
    ~ArrheniusRate();
    
    void set_Cf( const CoeffType Cf );
    void set_eta( const CoeffType eta );
    void set_Ea( const CoeffType Ea );

    void scale_Ea( const CoeffType scale );

    CoeffType Cf() const;
    CoeffType eta() const;
    CoeffType Ea() const;

    //! \return the rate evaluated at \p T.
    template <typename StateType>
    StateType operator()(const StateType T) const;

    //! \return the derivative with respect to temperature evaluated at \p T.
    template <typename StateType>
    StateType derivative( const StateType T ) const;

    //! Simultaneously evaluate the rate and its derivative at \p T.
    template <typename StateType>
    void rate_and_derivative(const StateType T, StateType& rate, StateType& drate_dT) const;

    //! Formatted print, by default to \p std::cout
    void print(std::ostream& os = std::cout) const;

    //! Formatted print.
    friend std::ostream& operator<<(std::ostream& os, const ArrheniusRate& rate)
    {
      rate.print(os);
      return os;
    }

  private:

    CoeffType _Cf;
    CoeffType _eta;
    CoeffType _Ea;
    
  };

  template<typename CoeffType>
  ArrheniusRate<CoeffType>::ArrheniusRate(const CoeffType Cf, const CoeffType eta, const CoeffType Ea)
    : _Cf(Cf),
      _eta(eta),
      _Ea(Ea)
  {
    return;
  }

  template<typename CoeffType>
  ArrheniusRate<CoeffType>::~ArrheniusRate()
  {
    return;
  }

  template<typename CoeffType>
  void ArrheniusRate<CoeffType>::print(std::ostream& os) const
  {
    os << _Cf;
    if (_eta != 0.) os << "*T^" << _eta;
    os << "*exp(-" << _Ea << "/T)";

    return;
  }

  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType>
  inline
  void ArrheniusRate<CoeffType>::set_Cf( const CoeffType Cf )
  {
    _Cf = Cf;
    return;
  }

  template<typename CoeffType>
  inline
  void ArrheniusRate<CoeffType>::set_eta( const CoeffType eta )
  {
    _eta = eta;
    return;
  }

  template<typename CoeffType>
  inline
  void ArrheniusRate<CoeffType>::set_Ea( const CoeffType Ea )
  {
    _Ea = Ea;
    return;
  }

  template<typename CoeffType>
  inline
  void ArrheniusRate<CoeffType>::scale_Ea( const CoeffType scale )
  {
    _Ea *= scale;
    return;
  }

  template<typename CoeffType>
  inline
  CoeffType ArrheniusRate<CoeffType>::Cf() const
  { return _Cf; }

  template<typename CoeffType>
  inline
  CoeffType ArrheniusRate<CoeffType>::eta() const
  { return _eta; }

  template<typename CoeffType>
  inline
  CoeffType ArrheniusRate<CoeffType>::Ea() const
  { return _Ea; }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType ArrheniusRate<CoeffType>::operator()(const StateType T) const
  {
    return _Cf* (std::pow(T,_eta)*std::exp(-_Ea/T));
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType ArrheniusRate<CoeffType>::derivative( const StateType T ) const
  {
    return (*this)(T)/T*(_eta + _Ea/T);
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  void ArrheniusRate<CoeffType>::rate_and_derivative( const StateType T,
						      StateType& rate,
						      StateType& drate_dT) const
  {
    rate     = (*this)(T);
    drate_dT = rate/T*(_eta + _Ea/T);
    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_ARRHENIUS_RATE_H
