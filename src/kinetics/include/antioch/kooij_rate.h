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

#ifndef ANTIOCH_KOOIJ_RATE_H
#define ANTIOCH_KOOIJ_RATE_H

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
  //! Kooij rate equation.
  /*!\class KooijRate
 *
   * The Kooij kinetics model is of the form:
   * \f[
   *   \alpha(T) =  A \left(\frac{T}{\mathrm{T_\text{ref}}}\right)^\beta \exp\left(-\frac{E_a}{\mathrm{R}T}\right)
   * \f]
   * thus
   * \f[
   *   \frac{\partial\alpha(T)}{\partial T} =  \alpha(T) \left(\frac{\beta}{T} + \frac{E_a}{\mathrm{R}T^2}\right)
   * \f]
   *
   * Internally, we use the reduced pre-exponential parameter \f$a = \frac{A}{\mathrm{T_\text{ref}}^\beta}\f$ and
   * the reduced activation energy \f$e = \frac{E_a}{\mathrm{R}}\f$
   */
  template<typename CoeffType=double>
  class KooijRate : public KineticsType<CoeffType>
  {

  private:
 
    CoeffType _raw_Cf;
    CoeffType _Cf;
    CoeffType _eta;
    CoeffType _raw_Ea;
    CoeffType _Ea;
    CoeffType _Tref;
    CoeffType _rscale;
  
  public:

    KooijRate (const CoeffType Cf=0., const CoeffType eta=0., const CoeffType Ea=0., const CoeffType Tref = 1., 
               const CoeffType rscale = Constants::R_universal<CoeffType>());
    ~KooijRate();
    
    void set_Cf(    const CoeffType Cf );
    void set_eta(   const CoeffType eta );
    void set_Ea(    const CoeffType Ea );
    void set_Tref(  const CoeffType Tref );
    void set_rscale(const CoeffType scale );

    //! set one parameter, characterized by enum
    void set_parameter(KineticsModel::Parameters parameter, const CoeffType & new_value);

    /*! reset the coeffs
     *
     * Two ways of modifying your rate:
     *   - you change totally the rate, thus you 
     *        require exactly five parameters, the order
     *        assumed is Cf, eta, Ea, Tref, rscale 
     *   - you just change the value, thus Tref and rscale are not
     *        modified. You require exactly three parameters,
     *        the order assumed is Cf, eta, Ea
     */
    template <typename VectorCoeffType>
    void reset_coefs(const VectorCoeffType & coefficients);

    CoeffType Cf()     const;
    CoeffType eta()    const;
    CoeffType Ea()     const;
    CoeffType Tref()   const;
    CoeffType rscale() const;

    //! \return the rate evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    rate(const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, _Cf* (ant_exp(_eta * ant_log(T) - _Ea/T)))

    //! \return the rate evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    operator()(const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, this->rate(T))

    //! \return the derivative with respect to temperature evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    derivative( const StateType& T ) const
    ANTIOCH_AUTOFUNC(StateType, (*this)(T)/T*(_eta + _Ea/T))

    //! Simultaneously evaluate the rate and its derivative at \p T.
    template <typename StateType>
    void rate_and_derivative(const StateType& T, StateType& rate, StateType& drate_dT) const;

// KineticsConditions overloads

    //! \return the rate evaluated at \p T.
    template <typename StateType, typename VectorStateType>
    ANTIOCH_AUTO(StateType) 
    rate(const KineticsConditions<StateType,VectorStateType>& T) const
    ANTIOCH_AUTOFUNC(StateType, _Cf * ant_exp(_eta * T.temp_cache().lnT - _Ea/T.T()))

    //! \return the rate evaluated at \p T.
    template <typename StateType, typename VectorStateType>
    ANTIOCH_AUTO(StateType) 
    operator()(const KineticsConditions<StateType,VectorStateType>& T) const
    ANTIOCH_AUTOFUNC(StateType, this->rate(T))

    //! \return the derivative with respect to temperature evaluated at \p T.
    template <typename StateType, typename VectorStateType>
    ANTIOCH_AUTO(StateType) 
    derivative( const KineticsConditions<StateType,VectorStateType>& T ) const
    ANTIOCH_AUTOFUNC(StateType, (*this)(T)/T.T()*(_eta + _Ea/T.T()))

    //! Simultaneously evaluate the rate and its derivative at \p T.
    template <typename StateType,typename VectorStateType>
    void rate_and_derivative(const KineticsConditions<StateType,VectorStateType>& T, StateType& rate, StateType& drate_dT) const;

    //! print equation
    const std::string numeric() const;

  private:

    void compute_cf();
   
  };

  template<typename CoeffType>
  KooijRate<CoeffType>::KooijRate(const CoeffType Cf, const CoeffType eta, const CoeffType Ea, const CoeffType Tref, const CoeffType rscale)
    : KineticsType<CoeffType>(KineticsModel::KOOIJ),
      _raw_Cf(Cf),
      _eta(eta),
      _raw_Ea(Ea),
      _Tref(Tref),
      _rscale(rscale)
  {

    _Ea = _raw_Ea / _rscale;
    this->compute_cf();

    return;
  }

  template<typename CoeffType>
  KooijRate<CoeffType>::~KooijRate()
  {
    return;
  }

  template<typename CoeffType>
  const std::string KooijRate<CoeffType>::numeric() const
  {
    std::stringstream os;
    os << _raw_Cf;
    if (_eta != 0.) os << "*(T/" << _Tref << ")^" << _eta;
    os << "*exp(-" << _raw_Ea << "/(R*T))";

    return os.str();
  }

  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType>
  inline
  void KooijRate<CoeffType>::set_Cf( const CoeffType Cf )
  {
    _raw_Cf = Cf;
    this->compute_cf();

    return;
  }

  template<typename CoeffType>
  inline
  void KooijRate<CoeffType>::set_Tref( const CoeffType Tref )
  {
    _Tref = Tref;
    this->compute_cf();

    return;
  }

  template<typename CoeffType>
  inline
  void KooijRate<CoeffType>::set_eta( const CoeffType eta )
  {
    _eta = eta;
    this->compute_cf();
    return;
  }

  template<typename CoeffType>
  inline
  void KooijRate<CoeffType>::set_Ea( const CoeffType Ea )
  {
    _raw_Ea = Ea;
    _Ea = _raw_Ea / _rscale;
    return;
  }

  template<typename CoeffType>
  inline
  void KooijRate<CoeffType>::set_rscale( const CoeffType rscale )
  {
    _rscale = rscale;
    _Ea = _raw_Ea / _rscale;
    return;
  }

  template<typename CoeffType>
  template <typename VectorCoeffType>
  inline
  void KooijRate<CoeffType>::reset_coefs(const VectorCoeffType & coefficients)
  {
       // 3 or 5
     antioch_assert_greater(coefficients.size(),2);
     antioch_assert_less(coefficients.size(),6);
     antioch_assert_not_equal_to(coefficients.size(),4);

     if(coefficients.size() == 5)
     {
        this->set_rscale(coefficients[4]);
        this->set_Tref(coefficients[3]);
     }
     this->set_Cf(coefficients[0]);
     this->set_eta(coefficients[1]);
     this->set_Ea(coefficients[2]);
  }

  template<typename CoeffType>
  inline
  void KooijRate<CoeffType>::set_parameter(KineticsModel::Parameters parameter, const CoeffType & new_value)
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
        this->set_Ea(new_value);
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
  CoeffType KooijRate<CoeffType>::Cf() const
  { return _Cf; }

  template<typename CoeffType>
  inline
  CoeffType KooijRate<CoeffType>::eta() const
  { return _eta; }

  template<typename CoeffType>
  inline
  CoeffType KooijRate<CoeffType>::Ea() const
  { return _Ea; }

  template<typename CoeffType>
  inline
  CoeffType KooijRate<CoeffType>::Tref() const
  { return _Tref; }

  template<typename CoeffType>
  inline
  CoeffType KooijRate<CoeffType>::rscale() const
  { return _rscale; }

  template<typename CoeffType>
  template<typename StateType>
  inline
  void KooijRate<CoeffType>::rate_and_derivative( const StateType& T,
                                                  StateType& rate,
                                                  StateType& drate_dT) const
  {
    rate     = (*this)(T);
    drate_dT = rate/T*(_eta + _Ea/T);

    return;
  }

  template<typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void KooijRate<CoeffType>::rate_and_derivative( const KineticsConditions<StateType,VectorStateType>& T,
                                                  StateType& rate,
                                                  StateType& drate_dT) const
  {
    rate     = (*this)(T);
    drate_dT = rate/T.T()*(_eta + _Ea/T.T());

    return;
  }

  template<typename CoeffType>
  inline
  void KooijRate<CoeffType>::compute_cf()
  {
    _Cf = _raw_Cf * ant_pow(KineticsModel::Tref<CoeffType>()/_Tref,_eta);
    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_KOOIJ_RATE_H
