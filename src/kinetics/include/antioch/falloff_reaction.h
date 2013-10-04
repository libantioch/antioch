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

#ifndef ANTIOCH_FALLOFF_REACTION_H
#define ANTIOCH_FALLOFF_REACTION_H

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/cmath_shims.h"
#include "antioch/reaction.h"
#include "antioch/lindemann_falloff.h"
#include "antioch/troe_falloff.h"

//C++
#include <string>
#include <vector>
#include <iostream>

namespace Antioch
{
  /*!\class FalloffReaction
  * Base class for falloff processes
  *
  * This class encapsulates a falloff reaction.  It performs the common operations.
  * A falloff rate constant is defined by the equation
  * \f[
  *     k(T,[M]) = \frac{[M] k_0(T)}{1+[M]\frac{k_0(T)}{k_\infty(T)}}\times F
  * \f]
  * with 
  * \f[
  *  \begin{split}
  *     k_0(T)      & = \lim_{[M] \rightarrow 0} k(T,[M])\\ 
  *     k_\infty(T) & = \lim_{[M] \rightarrow \infty} k(T,[M])
  *  \end{split}
  * \f]
  * \f$k_0\f$ and \f$k_\infty\f$ being respectively the low and high pressure
  * rate constant limits, considered elementary (as pressure in these conditions is constant).
  * Thus 
  * \f[
  *   \begin{split}
  *     k_0(T)      & = \alpha_0(T) \\
  *     k_\infty(T) & = \alpha_\infty(T)
  *   \end{split}
  * \f] 
  * with \f$\alpha_{i,\,i=0,\infty}(T)\f$ any
  * kinetics model (see base class KineticsType), and \f$[M]\f$
  * the mixture concentration (or pressure, it's equivalent, \f$[M] = \frac{P}{\mathrm{R}T}\f$
  * in the ideal gas model).  All reactions are assumed to be reversible, the kinetics models are
  * assumed to be the same.
  *
  * We have:
  * \f[
  *     \begin{split}
  *       \frac{\partial k(T,[M])}{\partial T} & = k(T,[M]) F \left(
  *                                                                 \frac{\partial k_0(T)}{\partial T}\frac{1}{k_0(T)} -
  *                                                                 \frac{\partial k_0(T)}{\partial T}\frac{1}{k_0(T) + \frac{k_\infty(T)}{[M]}} +
  *                                                                 \frac{\partial k_\infty(T)}{\partial T}
  *                                                                     \frac{k_0(T)}{k_\infty \left(k_0(T) + \frac{k_\infty(T)}{[M]}\right)}
  *                                                           \right) +
  *                                                 k(T,[M]) \frac{\partial F}{\partial T} \\[10pt]
  *       \frac{\partial k(T,[M])}{\partial c_i} & = F \frac{k(T,[M])}{[M] + [M]^2\frac{k_0(T)}{k_\infty(T)}} +
  *                                                  k(T,[M]) \frac{\partial F}{\partial c_i}
  *     \end{split}
  * \f]
  *
  * By default, the falloff is LindemannFalloff and the kinetics model KooijRate.
  */
  template<typename CoeffType=double, typename FalloffType = LindemannFalloff<CoeffType> >
  class FalloffReaction: public Reaction<CoeffType>
  {
  public:

    //! Construct a single reaction mechanism.
    FalloffReaction( const unsigned int n_species,
                     const std::string &equation, 
                     const bool &reversible = true,
                     const ReactionType::ReactionType &falloffType = ReactionType::LINDEMANN_FALLOFF,
                     const KineticsModel::KineticsModel kin = KineticsModel::KOOIJ);
    
    virtual ~FalloffReaction();

    //!
    template <typename StateType, typename VectorStateType>
    StateType compute_forward_rate_coefficient( const VectorStateType& molar_densities,
                                                const StateType& T ) const;
    
    //!
    template <typename StateType, typename VectorStateType>
    void compute_forward_rate_coefficient_and_derivatives( const VectorStateType& molar_densities,
                                                           const StateType& T,  
                                                           StateType& kfwd,  
                                                           StateType& dkfwd_dT, 
                                                           VectorStateType& dkfkwd_dX) const;


    //! Return const reference to the falloff object
    const FalloffType &F() const;
    //! Return writeable reference to the falloff object
    FalloffType &F();

  protected:

    FalloffType _F;

  };

  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType, typename FalloffType>
  inline
  FalloffReaction<CoeffType,FalloffType>::FalloffReaction( const unsigned int n_species,
                                                           const std::string &equation,
                                                           const bool &reversible,
                                                           const ReactionType::ReactionType &falloffType, 
                                                           const KineticsModel::KineticsModel kin)
    :Reaction<CoeffType>(n_species,equation,reversible,falloffType,kin),
     _F(n_species)
     
  {}


  template<typename CoeffType, typename FalloffType>
  inline
  FalloffReaction<CoeffType,FalloffType>::~FalloffReaction()
  {
    return;
  }

  template<typename CoeffType, typename FalloffType>
  inline
  FalloffType &FalloffReaction<CoeffType,FalloffType>::F()
  {
    return _F;
  }

  template<typename CoeffType, typename FalloffType>
  inline
  const FalloffType &FalloffReaction<CoeffType,FalloffType>::F() const
  {
    return _F;
  }

  template<typename CoeffType, typename FalloffType>
  template<typename StateType, typename VectorStateType>
  inline
  StateType FalloffReaction<CoeffType,FalloffType>::compute_forward_rate_coefficient( const VectorStateType& molar_densities,
                                                                                      const StateType& T  ) const
  {
//falloff is k(T,[M]) = k0*[M]/(1 + [M]*k0/kinf) * F = k0 * ([M]^-1 + k0 * kinf^-1)^-1 * F    
    StateType M = Antioch::zero_clone(T);
    for(unsigned int i = 0; i < molar_densities.size(); i++)
    {
        M += molar_densities[i];
    }

    return (*this->_forward_rate[0])(T) / (ant_pow(M,-1) + (*this->_forward_rate[0])(T) /(*this->_forward_rate[1])(T)) * 
            _F(T,molar_densities,(*this->_forward_rate[0])(T),(*this->_forward_rate[1])(T));
  }

  template<typename CoeffType, typename FalloffType>
  template<typename StateType, typename VectorStateType>
  inline
  void FalloffReaction<CoeffType,FalloffType>::compute_forward_rate_coefficient_and_derivatives( const VectorStateType &molar_densities,
                                                                                                 const StateType& T,
                                                                                                 StateType& kfwd, 
                                                                                                 StateType& dkfwd_dT,
                                                                                                 VectorStateType& dkfwd_dX) const 
  {
    //variables, k0,kinf and derivatives
   StateType k0 = Antioch::zero_clone(T);
    StateType dk0_dT = Antioch::zero_clone(T);
    StateType kinf = Antioch::zero_clone(T);
    StateType dkinf_dT = Antioch::zero_clone(T);

    this->_forward_rate[0]->compute_rate_and_derivative(T,k0,dk0_dT);
    this->_forward_rate[1]->compute_rate_and_derivative(T,kinf,dkinf_dT);

    StateType M = Antioch::zero_clone(T);
    for(unsigned int i = 0; i < molar_densities.size(); i++)
    {
        M += molar_densities[i];
    }

    //F
    StateType f = Antioch::zero_clone(T);
    StateType df_dT = Antioch::zero_clone(T);
    VectorStateType df_dX = Antioch::zero_clone(molar_densities);
    _F.F_and_derivatives(T,molar_densities,k0,dk0_dT,kinf,dkinf_dT,f,df_dT,df_dX);

// k(T,[M]) = k0*[M]/(1 + [M]*k0/kinf) * F = k0 * ([M]^-1 + k0 * kinf^-1)^-1 * F    
    kfwd = k0 / (ant_pow(M,-1) + k0/kinf); //temp variable here for calculations dk_d{T,X}

//dk_dT = F * dkfwd_dT + kfwd * dF_dT
//      = F * kfwd * (dk0_dT/k0 - dk0_dT/(kinf/[M] + k0) + k0 * dkinf_dT/(kinf * (kinf/[M] + k0) ) )
//      + dF_dT * kfwd
    dkfwd_dT = f * kfwd * (dk0_dT/k0 - dk0_dT/(kinf/M + k0) + dkinf_dT * k0/(kinf * (kinf/M + k0)))
             + df_dT * kfwd;

    dkfwd_dX.resize(this->n_species(), kfwd);
//dkfwd_dX = F * dkfwd_dX + kfwd * dF_dX
//         = F * kfwd / ([M] +  [M]^2 k0/kinf) + kfwd * dF_dX
    for(unsigned int ic = 0; ic < this->n_species(); ic++)
      {
        dkfwd_dX[ic] = f * kfwd / (M + ant_pow(M,2) * k0/kinf) + df_dX[ic] * kfwd;
      }

    kfwd *= f; //finalize

    return;
  }
  
} // namespace Antioch

#endif // ANTIOCH_FALLOFF_REACTION_H
