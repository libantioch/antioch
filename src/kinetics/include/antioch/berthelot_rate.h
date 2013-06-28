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

#ifndef ANTIOCH_BERTHELOT_RATE_H
#define ANTIOCH_BERTHELOT_RATE_H

//Antioch
#include "antioch/kinetics_type.h"

// C++
#include <cmath>
#include <iostream>
#include <sstream>

namespace Antioch
{
  //! Berthelot rate equation.
  /*!
   * Berthelot rate equation.  Computes rates of the form
   * \f$ C_f\times \exp(D*T) \f$. 
   */
  template<typename CoeffType=double>
  class BerthelotRate:public KineticsType<CoeffType>
  {
  
  public:

    BerthelotRate (const CoeffType Cf=0., const CoeffType D=0.);
    ~BerthelotRate();
    
    void set_Cf( const CoeffType Cf );
    void set_D( const CoeffType D );

    void scale_D( const CoeffType scale );

    CoeffType Cf() const;
    CoeffType D() const;

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
    CoeffType _D;
    
  };

  template<typename CoeffType>
  BerthelotRate<CoeffType>::BerthelotRate(const CoeffType Cf, const CoeffType D)
    : KineticsType<CoeffType>(KineticsModel::BERTHELOT),
      _Cf(Cf),
      _D(D)
  {
    return;
  }

  template<typename CoeffType>
  BerthelotRate<CoeffType>::~BerthelotRate()
  {
    return;
  }

  template<typename CoeffType>
  const std::string BerthelotRate<CoeffType>::numeric() const
  {
    std::stringstream os;
    os << _Cf;
    os << "*exp(" << _D << "*T)";

    return os.str();
  }

  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType>
  inline
  void BerthelotRate<CoeffType>::set_Cf( const CoeffType Cf )
  {
    _Cf = Cf;
    return;
  }

  template<typename CoeffType>
  inline
  void BerthelotRate<CoeffType>::set_D( const CoeffType D )
  {
    _D = D;
    return;
  }

  template<typename CoeffType>
  inline
  void BerthelotRate<CoeffType>::scale_D( const CoeffType scale )
  {
    _D *= scale;
    return;
  }

  template<typename CoeffType>
  inline
  CoeffType BerthelotRate<CoeffType>::Cf() const
  { return _Cf; }

  template<typename CoeffType>
  inline
  CoeffType BerthelotRate<CoeffType>::D() const
  { return _D; }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType BerthelotRate<CoeffType>::operator()(const StateType& T) const
  {
    using std::exp;
    return _Cf* (exp(_D*T));
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType BerthelotRate<CoeffType>::rate(const StateType& T) const
  {
    using std::exp;
    return _Cf* (exp(_D*T));
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType BerthelotRate<CoeffType>::derivative( const StateType& T ) const
  {
    return (*this)(T)*_D;
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  void BerthelotRate<CoeffType>::rate_and_derivative( const StateType& T,
						      StateType& rate,
						      StateType& drate_dT) const
  {
    rate     = (*this)(T);
    drate_dT = rate*_D;
    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_BERTHELOT_RATE_H
