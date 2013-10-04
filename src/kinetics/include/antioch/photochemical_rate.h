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
#include "antioch/kinetics_enum.h"
#include "antioch/sigma_bin_converter.h"

//C++
#include <vector>

namespace Antioch{

  template<typename CoeffType, typename VectorCoeffType = std::vector<CoeffType> >
  class PhotochemicalRate:public KineticsType<CoeffType>
  {

     private:
       VectorCoeffType _cross_section;
       VectorCoeffType _lambda_grid;
       Coefftype _k;
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
       void calculate_rate_constant(const VectorStateType &hv_flux, const VectorStateType &hv_lambda);
       
       //! \return the rate
       CoeffType rate() const;

       //! \return the rate evaluated at \p T.
       CoeffType operator() const;

       //! \return the derivative with respect to temperature.
       CoeffType derivative() const;

       //! Simultaneously evaluate the rate and its derivative at \p T.
       template <typename StateType>
       void rate_and_derivative(StateType& rate, StateType& drate_dT) const;

       //! print equation
       const std::string numeric() const;

  };

  template<typename CoeffType, typename VectorCoeffType>
  inline
  PhotochemicalRate<CoeffType,VectorCoeffType>::PhotochemicalRate(const VectorCoeffType &cs, 
                                                                  const VectorCoeffType &lambda):
    KineticsType<CoeffType>(KineticsModel::PHOTOCHEM),
    _cross_section(cs),
    _lambda_grid(lambda),
    _k(-1.)
  {
    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  PhotochemicalRate<CoeffType,VectorCoeffType>::PhotochemicalRate():
    KineticsType<CoeffType>(KineticsModel::PHOTOCHEM),
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
                                                                             const VectorStateType &hv_lambda)
  {
//cross-section and lambda exists
     antioch_assert_greater(_cross_section.size(),0);
     antioch_assert_greater(_lambda_grid.size(),0);

     VectorCoeffType cross_section_on_hv;
      _converter.y_on_custom_old_grid(_lambda_grid,_cross_section,hv_lambda,cross_section_on_hv);
      Antioch::set_zero(_k);
      for(unsigned int ibin = 0; ibin < hv_lambda.size(); ibin++)
      {
          _k += cross_section_on_hv[ibin] * hv_flux[ibin];
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
  inline
  CoeffType PhotochemicalRate<CoeffType,VectorCoeffType>::rate() const
  {
     return _k;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  CoeffType PhotochemicalRate<CoeffType,VectorCoeffType>::operator() const
  {
     this->rate();
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  CoeffType PhotochemicalRate<CoeffType,VectorCoeffType>::derivative() const
  {
     return 0.L;
  }

  template<typename CoeffType, typename VectorCoeffType>
  template <typename StateType>
  inline
  void rate_and_derivative(StateType& rate, StateType& drate_dT) const
  {
    rate = _k;
    Antioch::set_zero(drate_dT);
    return;
  }

} //end namespace Antioch

#endif
