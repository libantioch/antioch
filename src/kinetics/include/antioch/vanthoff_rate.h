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

#ifndef ANTIOCH_VAN_T_HOFF_RATE_H
#define ANTIOCH_VAN_T_HOFF_RATE_H

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
  //! Van't Hoff rate equation.
  /*!\class VantHoffRate
   * 
   * The Van't Hoff kinetics model is of the form:
   * \f[
   *   \alpha(T) =  A \left(\frac{T}{\mathrm{T_\text{ref}}}\right)^\beta \exp\left(-\frac{E_a}{\mathrm{R}T} + D T\right)
   * \f]
   * thus
   * \f[
   *   \frac{\partial\alpha(T)}{\partial T} =  \alpha(T)\left(D + \frac{\beta}{T} + \frac{E_a}{\mathrm{R}T^2}\right)
   * \f]
   *
   * Internally, we use the reduced pre-exponential parameter \f$a = \frac{A}{\mathrm{T_\text{ref}}^\beta}\f$ and
   * the reduced activation energy \f$e = \frac{E_a}{\mathrm{R}}\f$
   */
  template<typename CoeffType=double>
  class VantHoffRate:public KineticsType<CoeffType>
  {

  private:
 
    CoeffType _raw_Cf;
    CoeffType _Cf;
    CoeffType _eta;
    CoeffType _raw_Ea;
    CoeffType _Ea;
    CoeffType _D;
    CoeffType _Tref;
    CoeffType _rscale;
  
  public:

    VantHoffRate (const CoeffType Cf=0., const CoeffType eta=0., const CoeffType Ea=0.,
                  const CoeffType D=0., const CoeffType Tref = 1.,
                  const CoeffType rscale = Constants::R_universal<CoeffType>());

    ~VantHoffRate();
    
    void set_Cf(    const CoeffType Cf );
    void set_eta(   const CoeffType eta );
    //! set _raw_Ea, unit is known (cal.mol-1, J.mol-1, whatever), _Ea is computed
    void set_Ea(    const CoeffType Ea );
    //! set _Ea, unit is K, _raw_Ea is computed
    void reset_Ea(    const CoeffType Ea );
    void set_D(     const CoeffType D );
    void set_Tref(  const CoeffType Tref );
    void set_rscale(const CoeffType rscale );

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
     *        require exactly six parameters, the order
     *        assumed is Cf, eta, Ea, D, Tref, rscale 
     *   - you just change the value, thus Tref and rscale are not
     *        modified. You require exactly four parameters,
     *        the order assumed is Cf, eta, Ea, D
     */
    template <typename VectorCoeffType>
    void reset_coefs(const VectorCoeffType & coefficients);

    CoeffType Cf()     const;
    CoeffType eta()    const;
    CoeffType Ea()     const;
    CoeffType Ea_K()   const;
    CoeffType D()      const;
    CoeffType Tref()   const;
    CoeffType rscale() const;

// fork methods

    //! \return the rate evaluated at \p T.
    template <typename InputType>
    ANTIOCH_AUTO(typename return_type<InputType>::type) 
    operator()(const InputType& T) const
    ANTIOCH_AUTOFUNC(typename return_type<InputType>::type, this->rate(T))

    //! \return the derivative with respect to the parameter evaluated at \p T.
    template <typename InputType>
    ANTIOCH_AUTO(typename return_type<InputType>::type) 
    sensitivity( const InputType& T, KineticsModel::Parameters par) const;

//

    //! \return the rate evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    rate(const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, _Cf* (ant_pow(T,_eta)*ant_exp(-_Ea/T + _D*T)))

    //! \return the derivative with respect to temperature evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    derivative( const StateType& T ) const
    ANTIOCH_AUTOFUNC(StateType, (*this)(T)*(_D + _eta/T + _Ea/(T*T)))

    //! Simultaneously evaluate the rate and its derivative at \p T.
    template <typename StateType>
    void rate_and_derivative(const StateType& T, StateType& rate, StateType& drate_dT) const;

    //! \return the derivative with respect to the parameter A evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    sensitivity_A( const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, ant_pow(T/_Tref,_eta) * ant_exp(-_Ea/T + _D * T))

    //! \return the derivative with respect to the parameter beta evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    sensitivity_beta( const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, ant_log(T/_Tref) * this->rate(T))

    //! \return the derivative with respect to the parameter Ea evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    sensitivity_Ea( const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, - this->rate(T) /(_rscale * T) )

    //! \return the derivative with respect to the parameter D evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    sensitivity_D( const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, T * this->rate(T) )

    //! \return the derivative with respect to the parameter Tref evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    sensitivity_Tref( const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, - _eta /_Tref * this->rate(T) )

    //! \return the derivative with respect to the parameter Rscale evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    sensitivity_Rscale( const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, - this->rate(T) * _Ea /(_rscale * T) ) // Ea = _raw_Ea / rscale


// KineticsConditions overloads

    //! \return the rate evaluated at \p T.
    template <typename StateType, typename VectorStateType>
    ANTIOCH_AUTO(StateType) 
    rate(const KineticsConditions<StateType,VectorStateType>& T) const
    ANTIOCH_AUTOFUNC(StateType, _Cf * ant_exp(_eta * T.temp_cache().lnT - _Ea/T.T() + _D*T.T()))

    //! \return the derivative with respect to temperature evaluated at \p T.
    template <typename StateType, typename VectorStateType>
    ANTIOCH_AUTO(StateType) 
    derivative( const KineticsConditions<StateType,VectorStateType>& T ) const
    ANTIOCH_AUTOFUNC(StateType, (*this)(T)*(_D + _eta/T.T() + _Ea/(T.temp_cache().T2)))

    //! Simultaneously evaluate the rate and its derivative at \p T.
    template <typename StateType, typename VectorStateType>
    void rate_and_derivative(const KineticsConditions<StateType,VectorStateType>& T, StateType& rate, StateType& drate_dT) const;

    //! \return the derivative with respect to the parameter A evaluated at \p T.
    template <typename StateType, typename VectorStateType>
    ANTIOCH_AUTO(StateType) 
    sensitivity_A( const KineticsConditions<StateType,VectorStateType>& cond) const
    ANTIOCH_AUTOFUNC(StateType, ant_pow(cond.T()/_Tref,_eta) * ant_exp(_D * cond.T()))

    //! \return the derivative with respect to the parameter beta evaluated at \p T.
    template <typename StateType, typename VectorStateType>
    ANTIOCH_AUTO(StateType)
    sensitivity_beta(const KineticsConditions<StateType,VectorStateType> & cond) const
    ANTIOCH_AUTOFUNC(StateType,this->rate(cond) * (cond.temp_cache().lnT - ant_log(_Tref)) )

    //! \return the derivative with respect to parameter Ea evaluated at \p T.
    template <typename StateType, typename VectorStateType>
    ANTIOCH_AUTO(StateType) 
    sensitivity_Ea( const KineticsConditions<StateType, VectorStateType> & cond ) const
    ANTIOCH_AUTOFUNC(StateType, - this->rate(cond) /(_rscale * cond.T()) )

    //! \return the derivative with respect to the parameter D evaluated at \p T.
    template <typename StateType, typename VectorStateType>
    ANTIOCH_AUTO(StateType) 
    sensitivity_D( const KineticsConditions<StateType,VectorStateType> & cond) const
    ANTIOCH_AUTOFUNC(StateType, this->rate(cond) * cond.T())

    //! \return the derivative with respect to the parameter Tref evaluated at \p T.
    template <typename StateType, typename VectorStateType>
    ANTIOCH_AUTO(StateType) 
    sensitivity_Tref( const KineticsConditions<StateType,VectorStateType> & cond) const
    ANTIOCH_AUTOFUNC(StateType, - _eta /_Tref * this->rate(cond) )

    //! \return the derivative with respect to the parameter rscale evaluated at \p T.
    template <typename StateType, typename VectorStateType>
    ANTIOCH_AUTO(StateType) 
    sensitivity_Rscale( const KineticsConditions<StateType,VectorStateType> & cond) const
    ANTIOCH_AUTOFUNC(StateType, - this->rate(cond) * _Ea /(_rscale * cond.T()) ) // Ea = _raw_Ea / rscale

    //! print equation
    const std::string numeric() const;

  private:

    void compute_cf();
   
  };

  template<typename CoeffType>
  VantHoffRate<CoeffType>::VantHoffRate(const CoeffType Cf, const CoeffType eta, const CoeffType Ea, const CoeffType D, const CoeffType Tref, const CoeffType rscale)
    : KineticsType<CoeffType>(KineticsModel::VANTHOFF),
      _raw_Cf(Cf),
      _eta(eta),
      _raw_Ea(Ea),
      _D(D),
      _Tref(Tref),
      _rscale(rscale)
  {
    _Ea = _raw_Ea / _rscale;
    this->compute_cf();

    return;
  }

  template<typename CoeffType>
  VantHoffRate<CoeffType>::~VantHoffRate()
  {
    return;
  }

  template<typename CoeffType>
  const std::string VantHoffRate<CoeffType>::numeric() const
  {
    std::stringstream os;
    os << _raw_Cf;
    if (_eta != 0.) os << "*(T/" << _Tref << ")^" << _eta;
    os << "*exp(-" << _raw_Ea << "/(R*T) + " << _D << "*T)";

    return os.str();
  }

  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType>
  inline
  void VantHoffRate<CoeffType>::set_Cf( const CoeffType Cf )
  {
    _raw_Cf = Cf;
    this->compute_cf();

    return;
  }

  template<typename CoeffType>
  inline
  void VantHoffRate<CoeffType>::set_Tref( const CoeffType Tref )
  {
    _Tref = Tref;
    this->compute_cf();

    return;
  }

  template<typename CoeffType>
  inline
  void VantHoffRate<CoeffType>::set_eta( const CoeffType eta )
  {
    _eta = eta;
    this->compute_cf();
    return;
  }

  template<typename CoeffType>
  inline
  void VantHoffRate<CoeffType>::set_Ea( const CoeffType Ea )
  {
    _raw_Ea = Ea;
    _Ea = _raw_Ea / _rscale;
    return;
  }

  template<typename CoeffType>
  inline
  void VantHoffRate<CoeffType>::reset_Ea( const CoeffType Ea )
  {
    _Ea = Ea;
    _raw_Ea = _Ea * _rscale;
    return;
  }

  template<typename CoeffType>
  inline
  void VantHoffRate<CoeffType>::set_rscale( const CoeffType rscale )
  {
    _rscale = rscale;
    _Ea = _raw_Ea / _rscale;
    return;
  }

  template<typename CoeffType>
  inline
  void VantHoffRate<CoeffType>::set_D( const CoeffType D )
  {
    _D = D;
    return;
  }

  template<typename CoeffType>
  template <typename VectorCoeffType>
  inline
  void VantHoffRate<CoeffType>::reset_coefs(const VectorCoeffType & coefficients)
  {
       // 4 or 6
     antioch_assert_greater(coefficients.size(),3);
     antioch_assert_less(coefficients.size(),7);
     antioch_assert_not_equal_to(coefficients.size(),5);

     if(coefficients.size() == 6)
     {
        this->set_rscale(coefficients[5]);
        this->set_Tref(coefficients[4]);
     }
     this->set_Cf(coefficients[0]);
     this->set_eta(coefficients[1]);
     this->set_Ea(coefficients[2]);
     this->set_D(coefficients[3]);
  }

  template<typename CoeffType>
  inline
  void VantHoffRate<CoeffType>::set_parameter(KineticsModel::Parameters parameter, CoeffType new_value)
  {
    switch(parameter)
    {
      case KineticsModel::Parameters::R_SCALE:
      {
        this->set_rscale(new_value);
      }
       break;
      case KineticsModel::Parameters::T_REF:
      {
        this->set_Tref(new_value);
      }
       break;
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
      case KineticsModel::Parameters::E:
      {
        this->reset_Ea(new_value);
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
  CoeffType VantHoffRate<CoeffType>::get_parameter(KineticsModel::Parameters parameter) const
  {
    switch(parameter)
    {
      case KineticsModel::Parameters::R_SCALE:
      {
         return this->rscale();
      }
       break;
      case KineticsModel::Parameters::T_REF:
      {
         return this->Tref();
      }
       break;
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
      case KineticsModel::Parameters::E:
      {
         return this->Ea();
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
  CoeffType VantHoffRate<CoeffType>::Cf() const
  { return _Cf; }

  template<typename CoeffType>
  inline
  CoeffType VantHoffRate<CoeffType>::eta() const
  { return _eta; }

  template<typename CoeffType>
  inline
  CoeffType VantHoffRate<CoeffType>::Ea() const
  { return _raw_Ea; }

  template<typename CoeffType>
  inline
  CoeffType VantHoffRate<CoeffType>::Ea_K() const
  { return _Ea; }

  template<typename CoeffType>
  inline
  CoeffType VantHoffRate<CoeffType>::D() const
  { return _D; }

  template<typename CoeffType>
  inline
  CoeffType VantHoffRate<CoeffType>::Tref() const
  { return _Tref; }

  template<typename CoeffType>
  inline
  CoeffType VantHoffRate<CoeffType>::rscale() const
  { return _rscale; }

  template<typename CoeffType>
  template<typename StateType>
  inline
  void VantHoffRate<CoeffType>::rate_and_derivative( const StateType& T,
                                                     StateType& rate,
                                                     StateType& drate_dT) const
  {
    rate     = (*this)(T);
    drate_dT = rate*(_D + _eta/T + _Ea/(T*T));
    return;
  }

  template<typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void VantHoffRate<CoeffType>::rate_and_derivative( const KineticsConditions<StateType,VectorStateType>& T,
                                                     StateType& rate,
                                                     StateType& drate_dT) const
  {
    rate     = (*this)(T);
    drate_dT = rate*(_D + _eta/T.T() + _Ea/(T.temp_cache().T2));
    return;
  }

  template<typename CoeffType>
  inline
  void VantHoffRate<CoeffType>::compute_cf()
  {
    _Cf = _raw_Cf * ant_pow(KineticsModel::Tref<CoeffType>()/_Tref,_eta);
    return;
  }

  template<typename CoeffType>
  template <typename InputType>
  inline
  ANTIOCH_AUTO(typename return_type<InputType>::type) 
  VantHoffRate<CoeffType>::sensitivity( const InputType& T, KineticsModel::Parameters par) const
  {
    switch(par)
    {
      case KineticsModel::Parameters::A:
      {
        return this->sensitivity_A(T);
      }
        break;
      case KineticsModel::Parameters::B:
      {
        return this->sensitivity_beta(T);
      }
        break;
      case KineticsModel::Parameters::E:
      {
        return this->sensitivity_Ea(T);
      }
        break;
      case KineticsModel::Parameters::D:
      {
        return this->sensitivity_D(T);
      }
        break;
      case KineticsModel::Parameters::T_REF:
      {
        return this->sensitivity_Tref(T);
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

#endif // ANTIOCH_VAN_T_HOFF_RATE_H
