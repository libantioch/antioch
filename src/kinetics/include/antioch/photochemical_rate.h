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

#ifndef ANTIOCH_PHOTOCHEMICAL_RATE_H
#define ANTIOCH_PHOTOCHEMICAL_RATE_H

//Antioch
#include "antioch/metaprogramming_decl.h"
#include "antioch/antioch_asserts.h"
#include "antioch/sigma_bin_converter.h"
#include "antioch/kinetics_type.h"
#include "antioch/particle_flux.h"

//C++
#include <vector>
#include <string>
#include <sstream>

namespace Antioch{
  /*!Photochemical rate
   *
   * \todo Need to find a place to store k once calculated,
   * and recalculate only if the photon flux has changed.
   * Might want also to re-put _cross_section_on_flux_grid,
   * as usually, changes will appear in the photon flux, not
   * the lambda grid.
   *
   */
  template<typename CoeffType, typename VectorCoeffType = std::vector<CoeffType> >
  class PhotochemicalRate:public KineticsType<CoeffType,VectorCoeffType>
  {

     private:
       VectorCoeffType _cross_section;
       VectorCoeffType _lambda_grid;
       SigmaBinConverter<VectorCoeffType> _converter;

     public:
       PhotochemicalRate(const VectorCoeffType &cs, const VectorCoeffType &lambda);
       PhotochemicalRate();
       ~PhotochemicalRate();


       VectorCoeffType cross_section() const;
       VectorCoeffType lambda_grid()   const;


       //!
       void set_cross_section(const VectorCoeffType &cs);

       //!
       void set_lambda_grid(const VectorCoeffType &l);

       //!
       void set_cross_section(CoeffType cs, int il);

       //!
       void set_lambda_grid(CoeffType l, int il);

       //! set one value of one parameter, characterized by enum and its index
       void set_parameter(KineticsModel::Parameters parameter, int l, CoeffType new_value);

       //! set one parameter, characterized by enum
       CoeffType get_parameter(KineticsModel::Parameters parameter, int l) const;

       /*! reset the coefficients
        *
        *  Contrary to the other rate models, no choice here in the VectorCoeffType type.
        *  The order assumed is lambda, cross-section
        */
       void reset_coefs(const VectorCoeffType & coefficients);

       //! calculate _k for a given photon flux
       template<typename VectorStateType>
       void calculate_rate_constant(const VectorStateType &hv_flux, const VectorStateType &hv_lambda, bool x_update = true);
       
       //! \return the rate
       template <typename VectorStateType>
       ANTIOCH_AUTO(typename value_type<VectorStateType>::type)
          operator()(const ParticleFlux<VectorStateType> & pf) const
       ANTIOCH_AUTOFUNC(typename value_type<VectorStateType>::type , this->rate(pf))

       //! \return the rate evaluated at the given photon spectrum.
       template <typename VectorStateType>
       typename value_type<VectorStateType>::type 
                rate(const ParticleFlux<VectorStateType>& pf) const;

       //! \return the derivative with respect to temperature.
       template <typename VectorStateType>
       ANTIOCH_AUTO(typename value_type<VectorStateType>::type)
       derivative( const ParticleFlux<VectorStateType>& /* ex */) const
       ANTIOCH_AUTOFUNC(typename value_type<VectorStateType>::type, 0)

       //! Simultaneously evaluate the rate and its derivative at \p T.
       template <typename StateType, typename VectorStateType>
       void rate_and_derivative(const ParticleFlux<VectorStateType>& pf, StateType& rate, StateType& drate_dT) const;

       //! print equation
       const std::string numeric() const;

  };

  template<typename CoeffType, typename VectorCoeffType>
  inline
  PhotochemicalRate<CoeffType,VectorCoeffType>::PhotochemicalRate(const VectorCoeffType &cs, 
                                                                  const VectorCoeffType &lambda):
    KineticsType<CoeffType,VectorCoeffType>(KineticsModel::PHOTOCHEM),
    _cross_section(cs),
    _lambda_grid(lambda)
  {
    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  PhotochemicalRate<CoeffType,VectorCoeffType>::PhotochemicalRate():
    KineticsType<CoeffType,VectorCoeffType>(KineticsModel::PHOTOCHEM)
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
  inline
  void PhotochemicalRate<CoeffType,VectorCoeffType>::set_cross_section(CoeffType cs, int il)
  {
     antioch_assert_less(il,_cross_section.size());

    _cross_section[il] = cs;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  void PhotochemicalRate<CoeffType,VectorCoeffType>::set_lambda_grid(CoeffType l, int il)
  {
     antioch_assert_less(il,_lambda_grid.size());

     _lambda_grid[il] = l;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  void PhotochemicalRate<CoeffType,VectorCoeffType>::reset_coefs(const VectorCoeffType & coefficients)
  {
      // we require the two vectors being the same size
      antioch_assert_equal_to(coefficients.size()%2,0);

      const unsigned int subsize(coefficients.size()/2);
      VectorCoeffType l(subsize);
      VectorCoeffType cs(subsize);
      for(unsigned int il = 0; il < subsize; il++)
      {
          l[il] = coefficients[il];
          cs[il] = coefficients[il + subsize];
      }
      this->set_lambda_grid(l);
      this->set_cross_section(cs);
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  void PhotochemicalRate<CoeffType,VectorCoeffType>::set_parameter(KineticsModel::Parameters parameter, int l, CoeffType new_value)
  {

    switch(parameter)
    {
      case KineticsModel::Parameters::LAMBDA:
      {
         this->set_lambda_grid(new_value,l);
      }
       break;
      case KineticsModel::Parameters::SIGMA:
      {
         this->set_cross_section(new_value,l);
      }
        break;
      default:
      {
        antioch_error();
      }
       break;
    }
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  CoeffType PhotochemicalRate<CoeffType,VectorCoeffType>::get_parameter(KineticsModel::Parameters parameter,int l) const
  {

   antioch_assert_less(l,_cross_section.size());

    switch(parameter)
    {
      case KineticsModel::Parameters::LAMBDA:
      {
         return this->lambda_grid()[l];
      }
       break;
      case KineticsModel::Parameters::SIGMA:
      {
         return this->cross_section()[l];
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

  template<typename CoeffType, typename VectorCoeffType>
  template<typename VectorStateType>
  inline
  typename value_type<VectorStateType>::type 
        PhotochemicalRate<CoeffType,VectorCoeffType>::rate(const ParticleFlux<VectorStateType> & pf) const
  {
     const VectorStateType &hv_flux =  pf.flux();
     const VectorStateType &hv_lambda = pf.abscissa();

//cross-section and lambda exists
     antioch_assert_greater(_cross_section.size(),0);
     antioch_assert_greater(_lambda_grid.size(),0);

//needed grid to the right size
     VectorStateType cross_section_on_flux_grid(hv_lambda.size());

//put them on the right grid
      _converter.y_on_custom_grid(_lambda_grid,_cross_section,hv_lambda,cross_section_on_flux_grid);

//calculates
      typename value_type<VectorStateType>::type k;
      Antioch::set_zero(k);
      for(unsigned int ibin = 0; ibin < hv_lambda.size() - 1; ibin++)
      {
          k += cross_section_on_flux_grid[ibin] * hv_flux[ibin] * (hv_lambda[ibin+1] - hv_lambda[ibin]); //right stairs
      }
      return k;
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
  template <typename StateType, typename VectorStateType>
  inline
  void PhotochemicalRate<CoeffType,VectorCoeffType>::rate_and_derivative(const ParticleFlux<VectorStateType> &pf, StateType& rate, StateType& drate_dT) const
  {
    Antioch::init_clone(rate,this->rate(pf));
    Antioch::set_zero(drate_dT);
    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  VectorCoeffType PhotochemicalRate<CoeffType,VectorCoeffType>::cross_section() const
  {
     return _cross_section;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  VectorCoeffType PhotochemicalRate<CoeffType,VectorCoeffType>::lambda_grid() const
  {
     return _lambda_grid;
  }

} //end namespace Antioch

#endif
