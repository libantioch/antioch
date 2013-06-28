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

#ifndef ANTIOCH_THREEBODY_REACTION_H
#define ANTIOCH_THREEBODY_REACTION_H

// Antioch
#include "antioch/reaction.h"

//C++
#include <string>
#include <vector>
#include <iostream>

namespace Antioch
{
  //!A single reaction mechanism. 
  /*!
    This class encapsulates a three-body reaction process. A three-body process
    rate constant is defined by the equation
    \f[k(T,[M]) = \alpha(T)\times \sum_i\epsilon_iC_i\f]
    with \f$\alpha(T)\f$ being a kinetics model (see base class Reaction), \f$[M]\f$
    the mixture concentration (or pressure, it's equivalent, \f$[M] = \frac{P}{\mathrm{R}T}\f$
    in ideal gas model) and \f$C_i\f$ the concentration of species \f$i\f$.  All reactions are assumed to be reversible. 
    By default, the Kooij kinetics model is used.
  */
  template <typename CoeffType=double>
  class ThreeBodyReaction: public Reaction<CoeffType>
  {
  public:

    //! Construct a single reaction mechanism.
    ThreeBodyReaction( const unsigned int n_species, 
                       const std::string &equation,
                       const KineticsModel::KineticsModel kin = KineticsModel::KOOIJ);
    
    ~ThreeBodyReaction();

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
                                                           VectorStateType& dkfwd_dX) const;

  private:
    
  };

  /* ------------------------- Inline Functions -------------------------*/
  template <typename CoeffType>
  inline
  ThreeBodyReaction<CoeffType>::ThreeBodyReaction( const unsigned int n_species,
                                                   const std::string &equation ,
                                                   const KineticsModel::KineticsModel kin) 
    :Reaction<CoeffType>(n_species,equation,ReactionType::THREE_BODY,kin)
  {
    return;
  }

  template <typename CoeffType>
  inline
  ThreeBodyReaction<CoeffType>::~ThreeBodyReaction()
  {
    return;
  }



  template <typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  StateType ThreeBodyReaction<CoeffType>::compute_forward_rate_coefficient( const VectorStateType& molar_densities,
                                                                            const StateType& T  ) const
  {
    //k(T,[M]) = (sum eff_i * C_i) * ...
    StateType kfwd = (this->efficiency(0) * molar_densities[0] );

    for (unsigned int s=1; s<this->n_species(); s++)
      {        
        kfwd += ( this->efficiency(s) * molar_densities[s] );
      }

    //... alpha(T)
    return (kfwd * (*this->_forward_rate[0])(T));
  }


  template <typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void ThreeBodyReaction<CoeffType>::compute_forward_rate_coefficient_and_derivatives( const VectorStateType &molar_densities,
                                                                                       const StateType& T,  
                                                                                       StateType& kfwd, 
                                                                                       StateType& dkfwd_dT,
                                                                                       VectorStateType& dkfwd_dX) const
  {
    antioch_assert_equal_to(dkfwd_dX.size(),this->n_species());

    //dk_dT = dalpha_dT * [sum_s (eps_s * X_s)]
    //dk_dCi = alpha(T) * eps_i
    this->_forward_rate[0]->compute_rate_and_derivative(T,kfwd,dkfwd_dT);

    dkfwd_dX[0] = kfwd;
    StateType coef = (this->efficiency(0) * molar_densities[0]);

    for (unsigned int s=1; s<this->n_species(); s++)
      {        
        coef += ( this->efficiency(s) * molar_densities[s] );
        dkfwd_dX[s] = kfwd;
      }

    kfwd *= coef;
    dkfwd_dT *= coef;

    for (unsigned int s=0; s<this->n_species(); s++)
      {        
        dkfwd_dX[s] *= this->efficiency(s);
      }

    return;
  }
  
} // namespace Antioch

#endif // ANTIOCH_ELEMENTARY_REACTION_H
