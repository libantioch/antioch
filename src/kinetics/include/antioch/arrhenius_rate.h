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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef ANTIOCH_ARRHENIUS_RATE_H
#define ANTIOCH_ARRHENIUS_RATE_H

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/cmath_shims.h"
#include "antioch/kinetics_type.h"
#include "antioch/physical_constants.h"

// C++
#include <cmath>
#include <iostream>
#include <sstream>

namespace Antioch
{
  //! Arrhenius rate equation.
  /*!\class ArrheniusRate
 *
   * The Arrhenius kinetics model is of the form:
   * \f[
   *   \alpha(T) = A  \exp\left(-\frac{E_a}{\mathrm{R}T}\right) 
   * \f]
   * with \f$\mathrm{R}\f$ the ideal gas constant. We have:
   * \f[
   *   \frac{\partial\alpha(T)}{\partial T} = \alpha(T) \frac{E_a}{\mathrm{R}T^2}
   * \f]
   *
   * Internally, we use the reduced activation energy \f$e = \frac{E_a}{\mathrm{R}}\f$
   */
  template<typename CoeffType=double>
  class ArrheniusRate: public KineticsType<CoeffType>
  {
  
  // We declare private members early for use with decltype
  private:

    CoeffType _Cf;
    CoeffType _raw_Ea;
    CoeffType _Ea;
    CoeffType _rscale;
    
  public:

    ArrheniusRate (const CoeffType Cf=0., const CoeffType Ea=0., const CoeffType rscale = Constants::R_universal<CoeffType>());
    ~ArrheniusRate();
    
    void set_Cf(     const CoeffType Cf );
    //! set Ea, rescale the value, unit is known
    void set_Ea(     const CoeffType Ea );
    //! set Ea, no rescaling, unit is K
    void reset_Ea(     const CoeffType Ea );
    void set_rscale( const CoeffType rscale );

    //! set one parameter, characterized by enum
    //
    // Beware of the Ea parameter, it \e must be in Kelvin
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
     *        require exactly three parameters, the order
     *        assumed is Cf, Ea, rscale 
     *   - you just change the value, thus rscale is not
     *        modified. You require exactly two parameters,
     *        the order assumed is Cf, Ea
     */
    template <typename VectorCoeffType>
    void reset_coefs(const VectorCoeffType & coefficients);

    CoeffType Cf()     const;
    CoeffType Ea()     const;
    CoeffType Ea_K()   const;
    CoeffType rscale() const;

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
    ANTIOCH_AUTOFUNC(StateType, _Cf* (ant_exp(-_Ea/T)))

    //! \return the derivative with respect to temperature evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    derivative( const StateType& T ) const
    ANTIOCH_AUTOFUNC(StateType, (*this)(T)*(_Ea/(T*T)))

    //! \return the derivative with respect to the parameter A evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    sensitivity_A( const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, ant_exp(-_Ea / T))

    //! \return the derivative with respect to the parameter beta evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    sensitivity_Ea( const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, - this->rate(T) /(_rscale * T) )

    //! \return the derivative with respect to the parameter rscale evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    sensitivity_Rscale( const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, - this->rate(T) * _Ea /(_rscale * T) ) // Ea = _raw_Ea / rscale

    //! \return the derivative with respect to the parameter A evaluated at \p T.
    template <typename StateType, typename VectorStateType>
    ANTIOCH_AUTO(StateType) 
    sensitivity_A( const KineticsConditions<StateType, VectorStateType> & cond) const
    ANTIOCH_AUTOFUNC(StateType, ant_exp(-_Ea / cond.T()))

    //! \return the derivative with respect to the parameter beta evaluated at \p T.
    template <typename StateType, typename VectorStateType>
    ANTIOCH_AUTO(StateType) 
    sensitivity_Ea( const KineticsConditions<StateType, VectorStateType>& cond) const
    ANTIOCH_AUTOFUNC(StateType, - this->rate(cond.T()) /(_rscale * cond.T()) )

    //! \return the derivative with respect to the parameter rscale evaluated at \p T.
    template <typename StateType, typename VectorStateType>
    ANTIOCH_AUTO(StateType) 
    sensitivity_Rscale( const KineticsConditions<StateType, VectorStateType> & cond) const
    ANTIOCH_AUTOFUNC(StateType, - this->rate(cond.T()) * _Ea /(_rscale * cond.T()) ) // Ea = _raw_Ea / rscale

    //! Simultaneously evaluate the rate and its derivative at \p T.
    template <typename StateType>
    void rate_and_derivative(const StateType& T, StateType& rate, StateType& drate_dT) const;

    //! print equation
    const std::string numeric() const;

  };

  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType>
  ArrheniusRate<CoeffType>::ArrheniusRate(const CoeffType Cf, const CoeffType Ea, const CoeffType rscale)
    : KineticsType<CoeffType>(KineticsModel::ARRHENIUS),
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
  void ArrheniusRate<CoeffType>::reset_Ea( const CoeffType Ea )
  {
    _Ea = Ea;
    _raw_Ea = _Ea * _rscale;
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
  template <typename VectorCoeffType>
  inline
  void ArrheniusRate<CoeffType>::reset_coefs(const VectorCoeffType & coefficients)
  {
     antioch_assert_greater(coefficients.size(),1);
     antioch_assert_less(coefficients.size(),4);
     if(coefficients.size() == 3)this->set_rscale(coefficients[2]);
     this->set_Cf(coefficients[0]);
     this->set_Ea(coefficients[1]);
  }

  template<typename CoeffType>
  inline
  void ArrheniusRate<CoeffType>::set_parameter(KineticsModel::Parameters parameter, CoeffType new_value)
  {
      switch(parameter)
      {
        case KineticsModel::Parameters::A:
        {
          this->set_Cf(new_value);
        }
          break;
        case KineticsModel::Parameters::E:
        {
         this->reset_Ea(new_value);
        }
          break;
        case KineticsModel::Parameters::R_SCALE:
        {
         this->set_rscale(new_value);
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
  CoeffType ArrheniusRate<CoeffType>::get_parameter(KineticsModel::Parameters parameter) const
  {
      switch(parameter)
      {
        case KineticsModel::Parameters::A:
        {
          return this->Cf();
        }
          break;
        case KineticsModel::Parameters::E:
        {
         return this->Ea();
        }
          break;
        case KineticsModel::Parameters::R_SCALE:
        {
         return this->rscale();
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
  CoeffType ArrheniusRate<CoeffType>::Cf() const
  { return _Cf; }

  template<typename CoeffType>
  inline
  CoeffType ArrheniusRate<CoeffType>::Ea() const
  { return _raw_Ea; }

  template<typename CoeffType>
  inline
  CoeffType ArrheniusRate<CoeffType>::Ea_K() const
  { return _Ea; }

  template<typename CoeffType>
  inline
  CoeffType ArrheniusRate<CoeffType>::rscale() const
  { return _rscale; }

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

  template<typename CoeffType>
  template <typename InputType>
  inline
  ANTIOCH_AUTO(typename return_type<InputType>::type) 
  ArrheniusRate<CoeffType>::sensitivity( const InputType& T, KineticsModel::Parameters par ) const
  {
     switch(par)
     {
        case KineticsModel::Parameters::A:
        {
           return this->sensitivity_A(T);
        }
          break;
        case KineticsModel::Parameters::E:
        {
           return this->sensitivity_Ea(T);
        }
          break;
        case KineticsModel::Parameters::R_SCALE:
        {
           return this->sensitivity_Rscale(T);
        }
          break;
     }

    return typename return_type<InputType>::type(0);
  }

} // end namespace Antioch

#endif // ANTIOCH_ARRHENIUS_RATE_H
