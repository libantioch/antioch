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
#include "antioch/antioch_asserts.h"
#include "antioch/reaction.h"
#include "antioch/lindemann_falloff.h"

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
    By default, the falloff is Lindemann and the kinetics model Kooij.
  */
  template<typename CoeffType=double,typename FalloffType = LindemannFalloff<CoeffType> >
  class FalloffReaction: public Reaction<CoeffType>
  {
  public:

    //! Construct a single reaction mechanism.
    FalloffReaction( const unsigned int n_species,
                     const std::string &equation, 
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

    template <typename StateType, typename VectorStateType>
    StateType _Pr(const VectorStateType& molar_densities,
                                        const StateType& T ) const;

    template <typename StateType, typename VectorStateType>
    void _Pr_and_derivatives(const VectorStateType& molar_densities,
                                        const StateType& T ,
                                        const StateType& k0, 
                                        const StateType& kinf, 
                                        const StateType& dk0_dT, 
                                        const StateType& dkinf_dT, 
                                        StateType &Pr,
                                        StateType &dPr_dT,
                                        VectorStateType &dPr_dX) const;
    
  };

  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType, typename FalloffType>
  inline
  FalloffReaction<CoeffType,FalloffType>::FalloffReaction( const unsigned int n_species,
                                 const std::string &equation,
                                 const ReactionType::ReactionType &falloffType, 
                                 const KineticsModel::KineticsModel kin)
    :Reaction<CoeffType>(n_species,equation,falloffType,kin),
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
  StateType FalloffReaction<CoeffType,FalloffType>::_Pr( const VectorStateType& molar_densities,
                                                           const StateType& T  ) const
  {
//k(T,[M]) = [M] * alpha_0(T) /(1 + [M] * alpha_0(T)/alpha_inf(T)) * F
     StateType M = (molar_densities[0] );
     for (unsigned int s=1; s<this->n_species(); s++)
     {        
       M += molar_densities[s];
     }

//Pr = [M] k0 / kinf
     return (M * (*Reaction<CoeffType>::_forward_rate[0])(T) / (*Reaction<CoeffType>::_forward_rate[1])(T));
  }

  template<typename CoeffType, typename FalloffType>
  template<typename StateType, typename VectorStateType>
  inline
  void FalloffReaction<CoeffType,FalloffType>::_Pr_and_derivatives( const VectorStateType& molar_densities,
                                                           const StateType& T, 
                                                           const StateType& k0, 
                                                           const StateType& kinf, 
                                                           const StateType& dk0_dT, 
                                                           const StateType& dkinf_dT, 
                                                           StateType &Pr, 
                                                           StateType &dPr_dT,
                                                           VectorStateType &dPr_dX) const
  {
    antioch_assert_equal_to(dPr_dX.size(),this->n_species());
//k(T,[M]) = [M] * alpha_0(T) /(1 + [M] * alpha_0(T)/alpha_inf(T)) * F
    StateType M = (molar_densities[0] );
    for (unsigned int s=1; s<this->n_species(); s++)
    {        
      M += molar_densities[s];
    }
//Pr = [M] k0 / kinf
    Pr = M * (k0 / kinf);
    dPr_dT = M * dk0_dT / kinf - M * k0 * dkinf_dT / pow(kinf,2);
//dPr_dX = k0/kinf
    std::fill( dPr_dX.begin(),  dPr_dX.end(),  k0/kinf);
    return;   
  }

  template<typename CoeffType, typename FalloffType>
  template<typename StateType, typename VectorStateType>
  inline
  StateType FalloffReaction<CoeffType,FalloffType>::compute_forward_rate_coefficient( const VectorStateType& molar_densities,
                                                           const StateType& T  ) const
  {
    StateType Pr = _Pr(molar_densities,T);

    StateType one(1.);
    return (*this->_forward_rate[1])(T) * Pr /(one + Pr ) * _F(T,Pr);
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
    using std::pow;
//variables, k0,kinf and derivatives
    StateType k0 = Antioch::zero_clone(T);
    StateType dk0_dT = Antioch::zero_clone(T);
    StateType kinf = Antioch::zero_clone(T);
    StateType dkinf_dT = Antioch::zero_clone(T);
    this->_forward_rate[0]->compute_rate_and_derivative(T,k0,dk0_dT);
    this->_forward_rate[1]->compute_rate_and_derivative(T,kinf,dkinf_dT);

    StateType one = Antioch::zero_clone(T);
    one = 1.;
//Pr
    StateType Pr = Antioch::zero_clone(T);
    StateType dPr_dT = Antioch::zero_clone(T);
    VectorStateType dPr_dX = Antioch::zero_clone(molar_densities);
    _Pr_and_derivatives(molar_densities,T,k0,kinf,dk0_dT,dkinf_dT,Pr,dPr_dT,dPr_dX);

//F
    StateType f = Antioch::zero_clone(T);
    StateType df_dT = Antioch::zero_clone(T);
    VectorStateType df_dX = Antioch::zero_clone(molar_densities);
    _F.F_and_derivatives(T,Pr,dPr_dT,dPr_dX,f,df_dT,df_dX);

//k = kinf * Pr/(1 + Pr)
    kfwd = kinf * Pr / (one + Pr) * f;

//dk_dT = kfwd * [ dkinf_dT / kinf + dPr_dT /  Pr -  dPr_dT /(1 + Pr) + dF_dT / F]
    dkfwd_dT = kfwd * (dkinf_dT/kinf + dPr_dT / Pr - dPr_dT / (one + Pr) + df_dT/f);


    dkfwd_dX.resize(this->n_species(), kfwd);
//dkfwd_dX = kfwd * [dPr_dX/Pr - dPr_dX/(1+Pr) + dF_dX/F]
    for(unsigned int ic = 0; ic < this->n_species(); ic++)
    {
       dkfwd_dX[ic] = kfwd * (dPr_dX[ic]/(Pr*(one + Pr)) + df_dX[ic]/f) ;
    }

    return;
  }
  
} // namespace Antioch

#endif // ANTIOCH_FALLOFF_REACTION_H
