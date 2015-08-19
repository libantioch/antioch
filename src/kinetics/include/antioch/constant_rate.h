//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014 Paul T. Bauman, Benjamin S. Kirk, Sylvain Plessis,
//                    Roy H. Stonger
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

#ifndef ANTIOCH_CONSTANT_RATE_H
#define ANTIOCH_CONSTANT_RATE_H

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/kinetics_type.h"

// C++
#include <cmath>
#include <iostream>
#include <sstream>

namespace Antioch
{
  //! Constant rate equation.
  /*!\class ConstantRate
 *
   * This is a constant value for the kinetics model of the form:
   * \f[
   *  \alpha(T) = A 
   * \f]
   * Thus
   * \f[
   *  \frac{\partial \alpha(T)}{\partial T} = 0
   * \f]
   * 
   */
  template<typename CoeffType=double>
  class ConstantRate : public KineticsType<CoeffType>
  {

  private:
  
    CoeffType _Cf;
    
  public:

    ConstantRate (const CoeffType Cf=0.);
    ~ConstantRate();
    
    void set_Cf(  const CoeffType Cf );

    /*! reset the coeffs
     *
     * You require exactly one parameter
     */
    template <typename VectorCoeffType>
    void reset_coefs(const VectorCoeffType & coefficients);

    CoeffType Cf()   const;

    //! set one parameter, characterized by enum
    void set_parameter(KineticsModel::Parameters parameter, CoeffType new_value);

    //! get one parameter, characterized by enum
    CoeffType get_parameter(KineticsModel::Parameters parameter) const;

    //! for compatibility purpose with photochemistry (particle flux reactions)
    //
    // \todo, solve this
    template <typename VectorCoeffType>
    void set_parameter(KineticsModel::Parameters parameter, VectorCoeffType new_value){antioch_error();}

    //! \return the rate evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    rate(const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, constant_clone(T,_Cf))

    //! \return the rate evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    operator()(const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, this->rate(T))

    //! \return the derivative with respect to temperature evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    derivative( const StateType& T ) const
    ANTIOCH_AUTOFUNC(StateType, zero_clone(T))

    //! Simultaneously evaluate the rate and its derivative at \p T.
    template <typename StateType>
    void rate_and_derivative(const StateType& T, StateType& rate, StateType& drate_dT) const;

   // KineticsConditions overloads
    //! \return the rate evaluated at \p T.
    template <typename StateType, typename VectorStateType>
    ANTIOCH_AUTO(StateType) 
    rate(const KineticsConditions<StateType,VectorStateType>& cond) const
    ANTIOCH_AUTOFUNC(StateType, constant_clone(cond.T(),_Cf))

    //! \return the rate evaluated at \p T.
    template <typename StateType, typename VectorStateType>
    ANTIOCH_AUTO(StateType) 
    operator()(const KineticsConditions<StateType,VectorStateType>& cond) const
    ANTIOCH_AUTOFUNC(StateType, this->rate(cond))

    //! \return the derivative with respect to temperature evaluated at \p T.
    template <typename StateType, typename VectorStateType>
    ANTIOCH_AUTO(StateType) 
    derivative( const KineticsConditions<StateType,VectorStateType>& cond ) const
    ANTIOCH_AUTOFUNC(StateType, zero_clone(cond.T()))

    //! Simultaneously evaluate the rate and its derivative at \p T.
    template <typename StateType, typename VectorStateType>
    void rate_and_derivative(const KineticsConditions<StateType,VectorStateType>& cond, StateType& rate, StateType& drate_dT) const;

    //! print equation
    const std::string numeric() const;

  };

  template<typename CoeffType>
  ConstantRate<CoeffType>::ConstantRate(const CoeffType Cf)
    : KineticsType<CoeffType>(KineticsModel::CONSTANT),
      _Cf(Cf)
  {
    return;
  }

  template<typename CoeffType>
  ConstantRate<CoeffType>::~ConstantRate()
  {
    return;
  }

  template<typename CoeffType>
  const std::string ConstantRate<CoeffType>::numeric() const
  {
    std::stringstream os;
    os << _Cf;

    return os.str();
  }

  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType>
  inline
  void ConstantRate<CoeffType>::set_Cf( const CoeffType Cf )
  {
    _Cf = Cf;

    return;
  }

  template<typename CoeffType>
  inline
  void ConstantRate<CoeffType>::set_parameter(KineticsModel::Parameters parameter, CoeffType new_value)
  {
    antioch_assert_equal_to(parameter,KineticsModel::Parameters::A);

    this->set_Cf(new_value);

    return;
  }

  template<typename CoeffType>
  inline
  CoeffType ConstantRate<CoeffType>::get_parameter(KineticsModel::Parameters parameter) const
  {
    antioch_assert_equal_to(parameter,KineticsModel::Parameters::A);

    return this->Cf();
  }

  template<typename CoeffType>
  template <typename VectorCoeffType>
  inline
  void ConstantRate<CoeffType>::reset_coefs(const VectorCoeffType & coefficients)
  {
     antioch_assert_equal_to(coefficients.size(),1);

     this->set_Cf(coefficients[0]);
  }

  template<typename CoeffType>
  inline
  CoeffType ConstantRate<CoeffType>::Cf() const
  { return _Cf; }

  template<typename CoeffType>
  template<typename StateType>
  inline
  void ConstantRate<CoeffType>::rate_and_derivative( const StateType& /*T*/,
                                                          StateType& rate,
                                                          StateType& drate_dT) const
  {
    Antioch::constant_fill(rate, _Cf);
    Antioch::set_zero(drate_dT);
    return;
  }

  template<typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void ConstantRate<CoeffType>::rate_and_derivative(const KineticsConditions<StateType,VectorStateType>& /*cond*/, 
                                                    StateType& rate, 
                                                    StateType& drate_dT) const
  {
    Antioch::constant_fill(rate, _Cf);
    Antioch::set_zero(drate_dT);
    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_HERCOURT_ESSEN_RATE_H
