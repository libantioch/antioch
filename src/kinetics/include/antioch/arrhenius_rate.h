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
#include "antioch/physical_constants.h"

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

    ArrheniusRate (const CoeffType Cf=0., const CoeffType Ea=0., const CoeffType rscale = Constants::R_universal<CoeffType>()/1000.);
    ~ArrheniusRate();
    
    void set_Cf(     const CoeffType Cf );
    void set_Ea(     const CoeffType Ea );
    void set_rscale( const CoeffType rscale );

    CoeffType Cf()     const;
    CoeffType Ea()     const;
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

    CoeffType _Cf;
    CoeffType _raw_Ea;
    CoeffType _Ea;
    CoeffType _rscale;
    
  };

  template<typename CoeffType>
  ArrheniusRate<CoeffType>::ArrheniusRate(const CoeffType Cf, const CoeffType Ea, const CoeffType rscale)
    : KineticsType<CoeffType>(KinMod::ARRHENIUS),
      _Cf(Cf),
      _raw_Ea(Ea),
      _rscale(rscale)
  {
    _Ea = _raw_Ea / _rscale;
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
    os << "*exp(-" << _raw_Ea << "/(R*T))";

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
    _raw_Ea = Ea;
    _Ea = _raw_Ea / _rscale;
    return;
  }

  template<typename CoeffType>
  inline
  void ArrheniusRate<CoeffType>::set_rscale( const CoeffType rscale )
  {
    _rscale = rscale;
    _Ea = _raw_Ea / _rscale;
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
  inline
  CoeffType ArrheniusRate<CoeffType>::rscale() const
  { return _rscale; }

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
