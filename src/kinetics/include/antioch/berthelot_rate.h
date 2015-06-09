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

#ifndef ANTIOCH_BERTHELOT_RATE_H
#define ANTIOCH_BERTHELOT_RATE_H

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
  //! Berthelot rate equation.
  /*!\class BerthelotRate
 *
   * The Berthelot kinetics model is of the form:
   * \f[
   *    \alpha(T) =  A \exp\left(D T\right)
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
    
    void set_Cf( const CoeffType Cf );
    void set_D( const CoeffType D );

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
     * You require exactly two parameters, the order assumed is Cf, D
     */
    template <typename VectorCoeffType>
    void reset_coefs(const VectorCoeffType & coefficients);

    // \todo delete this method, why does it exist?
    void scale_D( const CoeffType scale );

    CoeffType Cf() const;
    CoeffType D() const;

    //! \return the rate evaluated at \p T.
    template <typename InputType>
    ANTIOCH_AUTO(typename return_type<InputType>::type) 
    operator()(const InputType& T) const
    ANTIOCH_AUTOFUNC(typename return_type<InputType>::type, this->rate(T))

    //! \return the derivative with respect to the parameter evaluated at \p T.
    template <typename InputType>
    ANTIOCH_AUTO(typename return_type<InputType>::type) 
    sensitivity( const InputType& T, KineticsModel::Parameters par ) const;

//
    //! \return the rate evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    rate(const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, _Cf* (ant_exp(_D*T)))

    //! \return the derivative with respect to temperature evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    derivative( const StateType& T ) const
    ANTIOCH_AUTOFUNC(StateType, (*this)(T)*_D)

    //! \return the derivative with respect to the parameter A evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    sensitivity_A( const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, ant_exp(_D * T))

    //! \return the derivative with respect to the parameter beta evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    sensitivity_D( const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, this->rate(T) * T)

    //! Simultaneously evaluate the rate and its derivative at \p T.
    template <typename StateType>
    void rate_and_derivative(const StateType& T, StateType& rate, StateType& drate_dT) const;

    //! \return the derivative with respect to the parameter A evaluated at \p T.
    template <typename StateType, typename VectorStateType>
    ANTIOCH_AUTO(StateType) 
    sensitivity_A( const KineticsConditions<StateType,VectorStateType> & cond) const
    ANTIOCH_AUTOFUNC(StateType, ant_exp(_D * cond.T()))

    //! \return the derivative with respect to the parameter beta evaluated at \p T.
    template <typename StateType, typename VectorStateType>
    ANTIOCH_AUTO(StateType) 
    sensitivity_D( const KineticsConditions<StateType,VectorStateType> & cond) const
    ANTIOCH_AUTOFUNC(StateType, this->rate(cond.T()) * cond.T())

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
  template <typename VectorCoeffType>
  inline
  void BerthelotRate<CoeffType>::reset_coefs(const VectorCoeffType & coefficients)
  {
     antioch_assert_equal_to(coefficients.size(),2);
     this->set_Cf(coefficients[0]);
     this->set_D(coefficients[1]);
  }

  template<typename CoeffType>
  inline
  void BerthelotRate<CoeffType>::set_parameter(KineticsModel::Parameters parameter, CoeffType new_value)
  {
     switch(parameter)
     {
       case KineticsModel::Parameters::A:
       {
         this->set_Cf(new_value);
       }
         break;
       case KineticsModel::Parameters::D:
       {
         this->set_D(new_value);
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
  CoeffType BerthelotRate<CoeffType>::get_parameter(KineticsModel::Parameters parameter) const
  {
     switch(parameter)
     {
       case KineticsModel::Parameters::A:
       {
         return this->Cf();
       }
         break;
       case KineticsModel::Parameters::D:
       {
         return this->D();
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

  template<typename CoeffType>
  template <typename InputType>
  inline
  ANTIOCH_AUTO(typename return_type<InputType>::type) 
  BerthelotRate<CoeffType>::sensitivity( const InputType& T, KineticsModel::Parameters par ) const
  {
     switch(par)
     {
        case KineticsModel::Parameters::A:
        {
           return this->sensitivity_A(T);
        }
          break;
        case KineticsModel::Parameters::D:
        {
           return this->sensitivity_D(T);
        }
          break;
     }

    return typename return_type<InputType>::type(0);
  }

} // end namespace Antioch

#endif // ANTIOCH_BERTHELOT_RATE_H
