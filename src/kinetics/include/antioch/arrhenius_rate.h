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
    void set_Ea(     const CoeffType Ea );
    void set_rscale( const CoeffType rscale );

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
    CoeffType rscale() const;

    //! \return the rate evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    rate(const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, _Cf* (ant_exp(-_Ea/T)))

    //! \return the rate evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    operator()(const StateType& T) const
    ANTIOCH_AUTOFUNC(StateType, this->rate(T))

    //! \return the derivative with respect to temperature evaluated at \p T.
    template <typename StateType>
    ANTIOCH_AUTO(StateType) 
    derivative( const StateType& T ) const
    ANTIOCH_AUTOFUNC(StateType, (*this)(T)*(_Ea/(T*T)))

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
