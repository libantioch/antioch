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

//Antioch
#include "antioch/kinetics_type.h"

// C++
#include <cmath>
#include <iostream>
#include <sstream>


namespace Antioch
{
  //! Arrhenius rate equation.
  /*!
   * Arrhenius rate equation.  Computes rates of the form
   * \f$ C_f\times \exp(-E_a/T) \f$. This class copied from
   * the \p FIN-S code and slightly reformatted for \p Antioch.
   */
  template<typename CoeffType=double>
  class ArrheniusRate: public KineticsType<CoeffType>
  {
  
  public:

    ArrheniusRate (const CoeffType Cf=0., const CoeffType Ea=0.);
    ~ArrheniusRate();
    
    void set_Cf( const CoeffType Cf );
    void set_Ea( const CoeffType Ea );

    void scale_Ea( const CoeffType scale );

    CoeffType Cf() const;
    CoeffType Ea() const;

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

    CoeffType _Cf;
    CoeffType _Ea;
    
  };

  template<typename CoeffType>
  ArrheniusRate<CoeffType>::ArrheniusRate(const CoeffType Cf, const CoeffType Ea)
    : KineticsType<CoeffType>(KinMod::ARRHENIUS),
      _Cf(Cf),
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
  const std::string ArrheniusRate<CoeffType>::numeric() const
  {
    std::stringstream os;
    os << _Cf;
    os << "*exp(-" << _Ea << "/T)";

    return os.str();
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
  CoeffType ArrheniusRate<CoeffType>::Ea() const
  { return _Ea; }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType ArrheniusRate<CoeffType>::operator()(const StateType& T) const
  {
    using std::exp;
    return _Cf* (exp(-_Ea/T));
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType ArrheniusRate<CoeffType>::rate(const StateType& T) const
  {
    using std::exp;
    return _Cf* (exp(-_Ea/T));
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType ArrheniusRate<CoeffType>::derivative( const StateType& T ) const
  {
    return (*this)(T)*(_Ea/(T*T));
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  void ArrheniusRate<CoeffType>::rate_and_derivative( const StateType& T,
						      StateType& rate,
						      StateType& drate_dT) const
  {
    rate     = (*this)(T);
    drate_dT = rate*_Ea/(T*T);
    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_ARRHENIUS_RATE_H
