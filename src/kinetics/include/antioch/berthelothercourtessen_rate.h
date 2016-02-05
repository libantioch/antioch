//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014-2016 Paul T. Bauman, Benjamin S. Kirk,
//                         Sylvain Plessis, Roy H. Stonger
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
#include "antioch/antioch_asserts.h"
#include "antioch/kinetics_type.h"
#include "antioch/cmath_shims.h"

// C++
#include <cmath>
#include <iostream>
#include <sstream>

namespace Antioch
{
  //! Berthelot Hercourt-Essen rate equation.
  /*!\class BerthelotHercourtEssenRate
 *
   * The Berthelot Hercourt Essen kinetics model is of the form:
   * \f[
   *    \alpha(T) =  A \left(\frac{T}{\mathrm{T_\text{ref}}}\right)^\beta  \exp\left(D T\right) 
   * \f]
   * thus
   * \f[
   *    \frac{\partial \alpha(T)}{\partial T} = \alpha(T) \left(\frac{\beta}{T} + D\right)
   * \f]
   *
   * Internally, we use the reduced temperature \f$t = \frac{T}{\mathrm{T_\text{ref}}}\f$.
   */
  template<typename CoeffType=double>
  class BerthelotHercourtEssenRate: public KineticsType<CoeffType>
  {

  private:
  
    CoeffType _raw_Cf;
    CoeffType _Cf;
    CoeffType _eta;
    CoeffType _D;
    CoeffType _Tref;
 
  public:

    BerthelotHercourtEssenRate (const CoeffType Cf=0., const CoeffType eta=0., const CoeffType D=0., const CoeffType Tref = 1.);
    ~BerthelotHercourtEssenRate();
    
    void set_Cf(  const CoeffType Cf );
    void set_eta( const CoeffType eta );
    void set_D(   const CoeffType D );
    void set_Tref(const CoeffType Tref );

    //! set one parameter, characterized by enum
    void set_parameter(KineticsModel::Parameters parameter, CoeffType new_value);

    //! get one parameter, characterized by enum
    CoeffType get_parameter(KineticsModel::Parameters parameter) const;

    //! for compatibility purpose with photochemistry (particle flux reactions)
    //
    // \todo, solve this
    template <typename VectorCoeffType>
    void set_parameter(KineticsModel::Parameters parameter, VectorCoeffType new_value){antioch_error();}

    /*! reset the coeffs
     *
     * Two ways of modifying your rate:
     *   - you change totally the rate, thus you 
     *        require exactly four parameters, the order
     *        assumed is Cf, eta, D, Tref
     *   - you just change the value, thus Tref is not
     *        modified. You require exactly three parameters,
     *        the order assumed is Cf, eta, D
     */
    template <typename VectorCoeffType>
    void reset_coefs(const VectorCoeffType & coefficients);

    CoeffType Cf()   const;
    CoeffType eta()  const;
    CoeffType D()    const;
    CoeffType Tref() const;

    //! \return the rate evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType)
    rate(const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, _Cf*(ant_pow(T,_eta)*ant_exp(_D*T)))

    //! \return the rate evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType)
    operator()(const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, this->rate(T))

    //! \return the derivative with respect to temperature evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType)
    derivative( const StateType& T ) const
    ANTIOCH_AUTOFUNC(StateType, (*this)(T)*(_eta/T + _D))

    //! Simultaneously evaluate the rate and its derivative at \p T.
    template <typename StateType>
    void rate_and_derivative(const StateType& T, StateType& rate, StateType& drate_dT) const;

//KineticsConditions overloads

    //! \return the rate evaluated at \p T.
    template <typename StateType, typename VectorStateType>
    ANTIOCH_AUTO(StateType)
    rate(const KineticsConditions<StateType,VectorStateType>& T) const
    ANTIOCH_AUTOFUNC(StateType, _Cf * ant_exp(_eta * T.temp_cache().lnT + _D*T.T()))

    //! \return the rate evaluated at \p T.
    template <typename StateType, typename VectorStateType>
    ANTIOCH_AUTO(StateType)
    operator()(const KineticsConditions<StateType,VectorStateType>& T) const
    ANTIOCH_AUTOFUNC(StateType, this->rate(T))

    //! \return the derivative with respect to temperature evaluated at \p T.
    template <typename StateType, typename VectorStateType>
    ANTIOCH_AUTO(StateType)
    derivative( const KineticsConditions<StateType,VectorStateType>& T ) const
    ANTIOCH_AUTOFUNC(StateType, (*this)(T)*(_eta/T.T() + _D))

    //! Simultaneously evaluate the rate and its derivative at \p T.
    template <typename StateType,typename VectorStateType>
    void rate_and_derivative(const KineticsConditions<StateType,VectorStateType>& T, StateType& rate, StateType& drate_dT) const;

    //! print equation
    const std::string numeric() const;

  private:

    void compute_cf();
  };

  template<typename CoeffType>
  BerthelotHercourtEssenRate<CoeffType>::BerthelotHercourtEssenRate(const CoeffType Cf, const CoeffType eta, const CoeffType D, const CoeffType Tref)
    : KineticsType<CoeffType>(KineticsModel::BHE),
      _raw_Cf(Cf),
      _eta(eta),
      _D(D),
      _Tref(Tref)
  {
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
    _raw_Cf = Cf;
    this->compute_cf();
    return;
  }

  template<typename CoeffType>
  inline
  void BerthelotHercourtEssenRate<CoeffType>::set_Tref( const CoeffType Tref )
  {
    _Tref = Tref;
    this->compute_cf();
    return;
  }

  template<typename CoeffType>
  inline
  void BerthelotHercourtEssenRate<CoeffType>::set_eta( const CoeffType eta )
  {
    _eta = eta;
    this->compute_cf();
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
  template <typename VectorCoeffType>
  inline
  void BerthelotHercourtEssenRate<CoeffType>::reset_coefs(const VectorCoeffType & coefficients)
  {
    antioch_assert_greater(coefficients.size(),2);
    antioch_assert_less(coefficients.size(),5);

    if(coefficients.size() == 4)this->set_Tref(coefficients[3]);
    this->set_Cf(coefficients[0]);
    this->set_eta(coefficients[1]);
    this->set_D(coefficients[2]);
  }

  template<typename CoeffType>
  inline
  void BerthelotHercourtEssenRate<CoeffType>::set_parameter(KineticsModel::Parameters parameter, CoeffType new_value)
  {
    switch(parameter)
    {
     case KineticsModel::Parameters::A:
     {
       this->set_Cf(new_value);
     }
      break;
     case KineticsModel::Parameters::B:
     {
       this->set_eta(new_value);
     }
      break;
     case KineticsModel::Parameters::D:
     {
       this->set_D(new_value);
     }
      break;
     case KineticsModel::Parameters::T_REF:
     {
       this->set_Tref(new_value);
     }
      break;
     default:
     {
       antioch_error();
     }
      break;
     }
  }

  template<typename CoeffType>
  inline
  CoeffType BerthelotHercourtEssenRate<CoeffType>::get_parameter(KineticsModel::Parameters parameter) const
  {
    switch(parameter)
    {
     case KineticsModel::Parameters::A:
     {
       return this->Cf();
     }
      break;
     case KineticsModel::Parameters::B:
     {
       return this->eta();
     }
      break;
     case KineticsModel::Parameters::D:
     {
       return this->D();
     }
      break;
     case KineticsModel::Parameters::T_REF:
     {
       return this->Tref();
     }
      break;
     default:
     {
       antioch_error();
     }
      break;
     }
     return 0;
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
  void BerthelotHercourtEssenRate<CoeffType>::rate_and_derivative( const StateType& T,
                                                                   StateType& rate,
                                                                   StateType& drate_dT) const
  {
    rate     = (*this)(T);
    drate_dT = rate*(_eta/T + _D);
    return;
  }

  template<typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void BerthelotHercourtEssenRate<CoeffType>::rate_and_derivative( const KineticsConditions<StateType,VectorStateType>& T,
                                                                   StateType& rate,
                                                                   StateType& drate_dT) const
  {
    rate     = (*this)(T);
    drate_dT = rate*(_eta/T.T() + _D);
    return;
  }

  template<typename CoeffType>
  inline
  void BerthelotHercourtEssenRate<CoeffType>::compute_cf()
  {
    _Cf = _raw_Cf * ant_pow(KineticsModel::Tref<CoeffType>()/_Tref,_eta);
    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_BERTHELOT_HERCOURT_ESSEN_RATE_H
