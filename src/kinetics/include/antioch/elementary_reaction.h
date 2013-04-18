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
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

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
    in ideal gas model).  All reactions are assumed to be reversible. 
  */
  class ElementaryReaction: public Reaction
  {
  public:

    //! Construct a single reaction mechanism.
    ElementaryReaction( const unsigned int n_species, const KinMod::KinMod kin, const std::string &equation );
    
    ~ElementaryReaction();

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

  private:
    
  };

  /* ------------------------- Inline Functions -------------------------*/
  inline
  ElementaryReaction::ElementaryReaction( const unsigned int n_species,
                                 const KinModel::KinModel kin,
				 const std::string &equation ) 
    :Reaction(n_species,kin,equation){}


  inline
  ElementaryReaction::~ElementaryReaction()
  {
    return;
  }



  template<typename StateType, typename VectorStateType>
  inline
  StateType ElementaryReaction::compute_forward_rate_coefficient( const VectorStateType& molar_densities,
							   const StateType& T  ) const
  {
//k(T,[M]) = alpha(T)
    return *_forward_rate[0](T);
  }

  template<typename StateType, typename VectorStateType>
  inline
  void ElementaryReaction::compute_forward_rate_coefficient_and_derivatives( const VectorStateType &molar_densities,
								      StateType& kfwd, 
								      StateType& dkfwd_dT,
								      VectorStateType& dkfwd_dY) const 
  {
//dk_dT = dalpha_dT(T)
     _forward_rate[0]->compute_rate_and_derivative(T,kfwd,dkfwd_dT);
//dk_dCi = 0.
    dkfwd_dY.resize(this->n_species(), kfwd);
    std::fill( dkfwd_dY.begin(),  dkfwd_dY.end(),  0.);    
    return;
  }
  
} // namespace Antioch

#endif // ANTIOCH_ELEMENTARY_REACTION_H
