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

#ifndef ANTIOCH_PHOTOCHEMICAL_RATE_H
#define ANTIOCH_PHOTOCHEMICAL_RATE_H

//Antioch
#include "antioch/metaprogramming.h"
#include "antioch/sigma_bin_converter.h"
#include "antioch/kinetics_type.h"

//C++
#include <vector>
#include <string>
#include <sstream>

namespace Antioch{

  template<typename CoeffType, typename VectorCoeffType = std::vector<CoeffType> >
  class PhotochemicalRate:public KineticsType<CoeffType,VectorCoeffType>
  {

     private:
       VectorCoeffType _cross_section;
       VectorCoeffType _lambda_grid;
       VectorCoeffType _cross_section_on_flux_grid;
       CoeffType _k;
       SigmaBinConverter<VectorCoeffType> _converter;

     public:
       PhotochemicalRate(const VectorCoeffType &cs, const VectorCoeffType &lambda);
       PhotochemicalRate();
       ~PhotochemicalRate();

       //!
       void set_cross_section(const VectorCoeffType &cs);

       //!
       void set_lambda_grid(const VectorCoeffType &l);

       //! calculate _k for a given photon flux
       template<typename VectorStateType>
       void calculate_rate_constant(const VectorStateType &hv_flux, const VectorStateType &hv_lambda, bool x_update = true);
       
       //! \return the rate
       template <typename StateType>
       ANTIOCH_AUTO(StateType) 
       operator()(const StateType& T) const
       ANTIOCH_AUTOFUNC(StateType, this->rate(T))

       //! \return the rate evaluated at \p T.
       template <typename StateType>
       ANTIOCH_AUTO(StateType) 
       rate(const StateType& T) const
       ANTIOCH_AUTOFUNC(StateType, Antioch::constant_clone(T, _k))

       //! \return the derivative with respect to temperature.
       template <typename StateType>
       ANTIOCH_AUTO(StateType) 
       derivative( const StateType& T ) const
       ANTIOCH_AUTOFUNC(StateType, Antioch::zero_clone(T))

       //! Simultaneously evaluate the rate and its derivative at \p T.
       template <typename StateType>
       void rate_and_derivative(const StateType& T, StateType& rate, StateType& drate_dT) const;

       //! print equation
       const std::string numeric() const;

  };

  template<typename CoeffType, typename VectorCoeffType>
  inline
  PhotochemicalRate<CoeffType,VectorCoeffType>::PhotochemicalRate(const VectorCoeffType &cs, 
                                                                  const VectorCoeffType &lambda):
    KineticsType<CoeffType,VectorCoeffType>(KineticsModel::PHOTOCHEM),
    _cross_section(cs),
    _lambda_grid(lambda),
    _k(-1.)
  {
    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  PhotochemicalRate<CoeffType,VectorCoeffType>::PhotochemicalRate():
    KineticsType<CoeffType,VectorCoeffType>(KineticsModel::PHOTOCHEM),
    _k(-1.)
  {
    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  PhotochemicalRate<CoeffType,VectorCoeffType>::~PhotochemicalRate()
  {
    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  void PhotochemicalRate<CoeffType,VectorCoeffType>::set_cross_section(const VectorCoeffType &cs)
  {
    _cross_section = cs;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  void PhotochemicalRate<CoeffType,VectorCoeffType>::set_lambda_grid(const VectorCoeffType &l)
  {
     _lambda_grid = l;
  }

  template<typename CoeffType, typename VectorCoeffType>
  template<typename VectorStateType>
  inline
  void PhotochemicalRate<CoeffType,VectorCoeffType>::calculate_rate_constant(const VectorStateType &hv_flux, 
                                                                             const VectorStateType &hv_lambda,
                                                                             bool x_update)
  {
//cross-section and lambda exists
     antioch_assert_greater(_cross_section.size(),0);
     antioch_assert_greater(_lambda_grid.size(),0);

      if(x_update || _cross_section_on_flux_grid.empty())
      {
        _cross_section_on_flux_grid.clear();
        _converter.y_on_custom_grid(_lambda_grid,_cross_section,hv_lambda,_cross_section_on_flux_grid);
      }
      Antioch::set_zero(_k);
      for(unsigned int ibin = 0; ibin < hv_lambda.size() - 1; ibin++)
      {
          _k += _cross_section_on_flux_grid[ibin] * hv_flux[ibin] * (hv_lambda[ibin+1] - hv_lambda[ibin]); //right stairs
      }
      return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  const std::string PhotochemicalRate<CoeffType,VectorCoeffType>::numeric() const
  {
    std::stringstream os;
    os << "int_0^infty sigma(lambda) * hv(lambda) * dlambda";

    return os.str();
  }


  template<typename CoeffType, typename VectorCoeffType>
  template <typename StateType>
  inline
  void PhotochemicalRate<CoeffType,VectorCoeffType>::rate_and_derivative(const StateType &T, StateType& rate, StateType& drate_dT) const
  {
    Antioch::constant_clone(rate,_k);
    Antioch::set_zero(drate_dT);
    return;
  }

} //end namespace Antioch

#endif
