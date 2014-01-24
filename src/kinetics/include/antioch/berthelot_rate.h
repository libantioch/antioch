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
#include "antioch/cmath_shims.h"

// C++
#include <cmath>
#include <iostream>
#include <sstream>

namespace Antioch
{
  //! Berthelot rate equation.
  /*!\class BerthelotRate
 *
   * The Berthelot kinetics model is of the form:
   * \f[
   *    \alpha(T) =  C_f \exp\left(D T\right)
   * \f] 
   * thus
   * \f[
   *    \frac{\partial\alpha(T)}{\partial T} = \alpha(T) D
   * \f] 
   */
  template<typename CoeffType=double>
  class BerthelotRate:public KineticsType<CoeffType>
  {

  private:

    CoeffType _Cf;
    CoeffType _D;
  
  public:

    BerthelotRate (const CoeffType Cf=0., const CoeffType D=0.);
    ~BerthelotRate();
    
    /** \brief set the pre-exponential factor in units {kmol,s,m} (e.g. order=1 : [s^-1]; order=2 : [m^3.kmol^-1.s^-1] )
     * 
     * \param Cf  pre-exponential factor in units {kmol,s,m} (e.g. order=1 : [s^-1]; order=2 : [m^3.kmol^-1.s^-1] )
     */ 
    void set_Cf( const CoeffType Cf );
    
    /** \brief set the exponential factor 
     * 
     * \param D  factor in [K^-1]
     */ 
    void set_D( const CoeffType D );

    void scale_D( const CoeffType scale );

    //! \return the pre-exponential factor in units {kmol,s,m} (e.g. order=1 : [s^-1]; order=2 : [m^3.kmol^-1.s^-1] )
    CoeffType Cf() const;
    
    //! \return the exponential factor in [K^-1]
    CoeffType D() const;

    //! \return the rate evaluated at \p T in units {kmol,s,m} (e.g. order=1 : [s^-1]; order=2 : [m^3.kmol^-1.s^-1] )..
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    rate(const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, _Cf* (ant_exp(_D*T)))

    //! \return the rate evaluated at \p T in units {kmol,s,m} (e.g. order=1 : [s^-1]; order=2 : [m^3.kmol^-1.s^-1] )..
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    operator()(const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, this->rate(T))

    //! \return the derivative with respect to temperature evaluated at \p T in units {kmol,s,m,K}.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    derivative( const StateType& T ) const
    ANTIOCH_AUTOFUNC(StateType, (*this)(T)*_D)

    //! Simultaneously evaluate the rate and its derivative at \p T in units {kmol,s,m,K}.
    template <typename StateType>
    void rate_and_derivative(const StateType& T, StateType& rate, StateType& drate_dT) const;

    //! print equation
    const std::string numeric() const;
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
