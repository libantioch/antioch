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

#ifndef ANTIOCH_ELEMENTARY_REACTION_H
#define ANTIOCH_ELEMENTARY_REACTION_H

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
    This class encapsulates an elementary reaction process. An elementary process
    rate constant is defined by the equation
    \f[k(T,[M]) = \alpha(T)\f]
    with \f$\alpha(T)\f$ being a kinetics model (see base class Reaction), and \f$[M]\f$
    the mixture concentration (or pressure, it's equivalent, \f$[M] = \frac{P}{\mathrm{R}T}\f$
    in ideal gas model).  All reactions are assumed to be reversible. By default, the kinetisc model
    is Kooij. 
  */
  template<typename CoeffType=double>
  class ElementaryReaction: public Reaction<CoeffType>
  {
  public:

    //! Construct a single reaction mechanism.
    ElementaryReaction( const unsigned int n_species, 
                        const std::string &equation,
                        const bool &reversible = true,
                        const KineticsModel::KineticsModel kin = KineticsModel::KOOIJ);
    
    ~ElementaryReaction();

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
  template<typename CoeffType>
  inline
  ElementaryReaction<CoeffType>::ElementaryReaction( const unsigned int n_species,
                                                     const std::string &equation,
                                                     const bool &reversible,
                                                     const KineticsModel::KineticsModel kin)
    :Reaction<CoeffType>(n_species,equation,reversible,ReactionType::ELEMENTARY,kin)
  {
    return;
  }


  template<typename CoeffType>
  inline
  ElementaryReaction<CoeffType>::~ElementaryReaction()
  {
    return;
  }



  template<typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  StateType ElementaryReaction<CoeffType>::compute_forward_rate_coefficient( const VectorStateType& molar_densities,
                                                                             const StateType& T) const
  {
    antioch_assert_equal_to(1, Reaction<CoeffType>::_forward_rate.size());

    //k(T,[M]) = alpha(T)
    return (*this->_forward_rate[0])(T);
  }

  template<typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void ElementaryReaction<CoeffType>::compute_forward_rate_coefficient_and_derivatives( const VectorStateType &molar_densities,
                                                                                        const StateType& T, 
                                                                                        StateType& kfwd, 
                                                                                        StateType& dkfwd_dT,
                                                                                        VectorStateType& dkfwd_dX) const 
  {
    //dk_dT = dalpha_dT(T)
    this->_forward_rate[0]->compute_rate_and_derivative(T,kfwd,dkfwd_dT);

    //dk_dCi = 0.
    antioch_assert_equal_to(dkfwd_dX.size(),this->n_species());
    Antioch::set_zero(dkfwd_dX);

    return;
  }
  
} // namespace Antioch

#endif // ANTIOCH_ELEMENTARY_REACTION_H
