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

#ifndef ANTIOCH_DUPLICATE_REACTION_H
#define ANTIOCH_DUPLICATE_REACTION_H

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
    This class encapsulates a duplicate reaction process. A duplicate process
    rate constant is defined by the equation
        \f[k(T,[M]) = \sum_i\alpha_i(T)\f]
    with \f$\alpha_i(T)\f$ being the \f$ith\f$ kinetics model (see base class Reaction), and \f$[M]\f$
    the mixture concentration (or pressure, it's equivalent, \f$[M] = \frac{P}{\mathrm{R}T}\f$
    in ideal gas model).  It is assumed that all the \f$\alpha_i\f$ are the same kinetisc model.
    All reactions are assumed to be reversible. By default, the kinetics model used is the Kooij equation.
  */
  template <typename CoeffType=double>
  class DuplicateReaction: public Reaction<CoeffType>
  {
  public:

    //! Construct a single reaction mechanism.
    DuplicateReaction( const unsigned int n_species, 
                       const std::string &equation,
                       const KinMod::KinMod kin = KinMod::KOOIJ);
    
    ~DuplicateReaction();

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
						   VectorStateType& dRbkwd_Y) const;

  private:
    
  };

  /* ------------------------- Inline Functions -------------------------*/
  template <typename CoeffType>
  inline
  DuplicateReaction<CoeffType>::DuplicateReaction( const unsigned int n_species,
				 const std::string &equation,
                                 const KinMod::KinMod kin)
    :Reaction<CoeffType>(n_species,equation,ReactionType::DUPLICATE,kin){}


  template <typename CoeffType>
  inline
  DuplicateReaction<CoeffType>::~DuplicateReaction()
  {
    return;
  }



  template <typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  StateType DuplicateReaction<CoeffType>::compute_forward_rate_coefficient( const VectorStateType& molar_densities,
							   const StateType& T  ) const
  {
    StateType kfwd = (*this->_forward_rate[0])(T);
    for(unsigned int ir = 1; ir < this->_forward_rate.size(); ir++)
    {
      kfwd += (*this->_forward_rate[ir])(T);
    }
    return kfwd;
  }

  template <typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void DuplicateReaction<CoeffType>::compute_forward_rate_coefficient_and_derivatives( const VectorStateType &molar_densities,
                                                                      const StateType& T, 
								      StateType& kfwd, 
								      StateType& dkfwd_dT,
						                      VectorStateType& dkfwd_dY) const
  {
//dk_dT = sum_p dalpha_p_dT
    StateType kfwd_tmp,dkfwd_dT_tmp;
     this->_forward_rate[0]->compute_rate_and_derivative(T,kfwd,dkfwd_dT);
    for(unsigned int ir = 1; ir < Reaction<CoeffType>::_forward_rate.size(); ir++)
    {
      this->_forward_rate[ir]->compute_rate_and_derivative(T,kfwd_tmp,dkfwd_dT_tmp);
      kfwd += kfwd_tmp;
      dkfwd_dT += dkfwd_dT_tmp;
    }

//dk_dCi = 0
    dkfwd_dY.resize(this->n_species(), kfwd);
    std::fill( dkfwd_dY.begin(),  dkfwd_dY.end(),  0.);    
    return;
  }
  
} // namespace Antioch

#endif // ANTIOCH_ELEMENTARY_REACTION_H
