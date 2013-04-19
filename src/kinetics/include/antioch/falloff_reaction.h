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
//
// $Id: elementary_reaction.h 38747 2013-04-17 23:26:39Z splessis $
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef ANTIOCH_FALLOFF_REACTION_H
#define ANTIOCH_FALLOFF_REACTION_H

// Antioch
#include "antioch/reaction.h"

//C++
#include <string>
#include <vector>
#include <iostream>

namespace Antioch
{
  //!Base class to falloff
  /*!
    This class encapsulates a falloff reaction.  It performs the common operations.
    A falloff rate constant is defined by the equation
        \f[k(T,[M]) = k_\infty(T)\frac{[M] k_0(T)}{1+[M]\frac{k_0(T)}{k_\infty(T)}}\times F\f]
    with \f$k_0(T) = k(T,[M] = [M]_0)\f$ and \f$k_\infty(T) = k(T,[M] = [M]_\infty)\f$; being respectively the low and high pressure
    rate constants, considered elementary (as pressure in these conditions is constant).  Thus \f$k_0(T) = \alpha_0(T)\f$ 
    and \f$k_\infty(T) = \alpha_\infty(T)\f$ with \f$\alpha_{i,\,i=0,\infty}(T)\$
    kinetics model (see base class Reaction), and \f$[M]\f$
    the mixture concentration (or pressure, it's equivalent, \f$[M] = \frac{P}{\mathrm{R}T}\f$
    in ideal gas model).  All reactions are assumed to be reversible, the kinetics models are
    assumed to be the same.
  */
  template<typename FalloffType>
  class FalloffReaction: public Reaction
  {
  public:

    //! Construct a single reaction mechanism.
    FalloffReaction( const unsigned int n_species, const KinMod::KinMod kin, const std::string &equation );
    
    virtual ~FalloffReaction();

    //!
    template <typename StateType, typename VectorStateType>
    StateType compute_forward_rate_coefficient( const VectorStateType& molar_densities,
					const StateType& kfwd, 
					const StateType& kbkwd ) const;
    
    //!
    template <typename StateType, typename VectorStateType>
    void compute_forward_rate_coefficient_and_derivatives( const VectorStateType& molar_densities,
						   StateType& kfwd,  
						   StateType& dkfwd_dT, 
						   VectorStateType& dkfkwd_dY) const;


    //! Return const reference to the falloff object
    const FalloffType &F() const;
    //! Return writeable reference to the falloff object
    FalloffType &F();

  protected:

    FalloffType _F;
    
  };

  /* ------------------------- Inline Functions -------------------------*/
  template<typename FalloffType>
  inline
  FalloffReaction<FalloffType>::FalloffReaction( const unsigned int n_species,
                                 const KinModel::KinModel kin,
				 const std::string &equation ) 
    :Reaction(n_species,kin,equation){}


  template<typename FalloffType>
  inline
  FalloffReaction<FalloffType>::~FalloffReaction()
  {
    return;
  }

  template<typename FalloffType>
  inline
  FalloffType &FalloffReaction<FalloffType>::F()
  {
     return _F;
  }

  template<typename FalloffType>
  inline
  const FalloffType &FalloffReaction<FalloffType>::F() const
  {
     return _F;
  }

  template<typename FalloffType>
  template<typename StateType, typename VectorStateType>
  inline
  StateType FalloffReaction::compute_forward_rate_coefficient( const VectorStateType& molar_densities,
							   const StateType& T  ) const
  {
//k(T,[M]) = [M] * alpha_0(T) /(1 + [M] * alpha_0(T)/alpha_inf(T)) * F
     StateType M = (molar_densities[0] );
     for (unsigned int s=1; s<this->n_species(); s++)
     {        
       M += molar_densities[s];
     }

//Pr = [M] k0 / kinf
     StateType Pr = M * (*_forward_rate[0])(T) / (*_forward_rate[1])(T)

     StateType kfwd = (*_forward_rate[1])(T) * Pr /(1. + Pr ) * _F.compute_F(T,Pr);
  }

  template<typename FalloffType>
  template<typename StateType, typename VectorStateType>
  inline
  void FalloffReaction::compute_forward_rate_coefficient_and_derivatives( const VectorStateType &molar_densities,
								      StateType& kfwd, 
								      StateType& dkfwd_dT,
								      VectorStateType& dkfwd_dY) const 
  {
    using std::pow;
//variables
    StateType k0,dk0_dT;
    _forward_rate[0]->rate_and_derivative(T,k0,dk0_dT);
    StateType kinf,dkinf_dT;
    _forward_rate[1]->rate_and_derivative(T,kinf,dkinf_dT);

    StateType M = (molar_densities[0] );
    for (unsigned int s=1; s<this->n_species(); s++)
    {        
      M += molar_densities[s];
    }
    StateType Pr = M * k0/kinf;
//k = kinf * Pr/(1 + Pr)
    kfwd = kinf * Pr / (1. + Pr);

    StateType dPr_dT = M * dk0_dT / kinf - M * k0 * dkinf_dT / pow(kinf,2);
//dk_dT = kfwd * [ dkinf_dT / kinf + dPr_dT /  Pr -  dPr_dT /(1 + Pr) ]
    dkfwd_dT = kfwd * (dkinf_dT/kinf + dPr_dT / Pr - dPr_dT / (1. + Pr));


    dkfwd_dY.resize(this->n_species(), kfwd);
//dM_dY = 1.
//dPr_dY = k0/kinf
//dkfwd_dY = kfwd * [dPr_dY/Pr - dPr_dY/(1+Pr)]
//dkfwd_dY = kfwd * [k0/kinf/Pr - k0/kinf/(1+Pr)]
    StateType dkfwd_dCi = kfwd * (k0/kinf/Pr - k0/kinf/(1. + Pr));
    std::fill( dkfwd_dY.begin(),  dkfwd_dY.end(),  dkfwd_dCi);    
    return;
  }
  
} // namespace Antioch

#endif // ANTIOCH_FALLOFF_REACTION_H
