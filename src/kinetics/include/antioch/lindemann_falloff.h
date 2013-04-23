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

#ifndef _ANTIOCH_LINDEMANN_FALLOFF_H
#define _ANTIOCH_LINDEMANN_FALLOFF_H

//Antioch

//C++

namespace Antioch{
/*!\class LindemannFalloff
 *
 * The Lindemann model is the simplest falloff
 * model:
 * \f[
 *     F = 1
 * \f]
 */
template <typename CoeffType = double>
class LindemannFalloff{
     public:
       LindemannFalloff();
       ~LindemannFalloff();

     template <typename StateType>
     StateType operator()(const StateType& T, const StateType &Pr) const;

     template <typename StateType, typename VectorStateType>
     StateType F_and_derivatives(const StateType& T, 
                           const StateType &Pr, 
                           const StateType &dPr_dT, 
                           const VectorStateType &dPr_dY,
                           StateType &F,
                           StateType &dF_dT,
                           VectorStateType &dF_dY) const;


};
  template<typename CoeffType>
  template <typename StateType>
  inline
  StateType LindemannFalloff<CoeffType>::operator()(const StateType &T, const StateType &Pr) const
  {
    return StateType(1.);
  }

  template <typename CoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  StateType LindemannFalloff<CoeffType>::F_and_derivatives(const StateType& T, 
                        const StateType &Pr, 
                        const StateType &dPr_dT, 
                        const VectorStateType &dPr_dY,
                        StateType &dF_dT,
                        StateType &F,
                        VectorStateType &dF_dY) const
  {
//all derived are 0
    dF_dT = 0.;
    dF_dY.resize(this->n_species(), 0.);
    std::fill( dF_dY.begin(),  dF_dY.end(),  0.);    
    return;
  }

  template<typename CoeffType>
  inline
  LindemannFalloff<CoeffType>::LindemannFalloff()
  {
    return;
  }

  template<typename CoeffType>
  inline
  LindemannFalloff<CoeffType>::~LindemannFalloff()
  {
    return;
  }
}

#endif
