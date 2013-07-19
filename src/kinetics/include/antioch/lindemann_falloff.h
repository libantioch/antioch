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
    LindemannFalloff(const unsigned int nspec);
    ~LindemannFalloff();

    template <typename StateType, typename VectorStateType>
    StateType operator()(const StateType &T,
                         const VectorStateType &molar_densities,
                         const StateType &k0, 
                         const StateType &kinf) const;

    template <typename StateType, typename VectorStateType>
    void F_and_derivatives(const StateType& T, 
                           const VectorStateType &molar_densities,
                           const StateType &k0, 
                           const StateType &dk0_dT, 
                           const StateType &kinf, 
                           const StateType &dkinf_dT, 
                           StateType &F,
                           StateType &dF_dT,
                           VectorStateType &dF_dX) const;

  private:
    unsigned int n_spec;

  };
  template<typename CoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  StateType LindemannFalloff<CoeffType>::operator()(const StateType &T,
                                                    const VectorStateType &molar_densities,
                                                    const StateType &k0, 
                                                    const StateType &kinf) const
  {
    return Antioch::constant_clone(T, 1);
  }

  template <typename CoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  void LindemannFalloff<CoeffType>::F_and_derivatives(const StateType& T, 
                                                      const VectorStateType &molar_densities,
                                                      const StateType &k0, 
                                                      const StateType &dk0_dT, 
                                                      const StateType &kinf, 
                                                      const StateType &dkinf_dT, 
                                                      StateType &F,
                                                      StateType &dF_dT,
                                                      VectorStateType &dF_dX) const
  {
    //all derived are 0
    Antioch::set_zero(dF_dT);
    antioch_assert_equal_to(dF_dX.size(),n_spec);
    Antioch::set_zero(dF_dX);
    // F = 1
    F = Antioch::constant_clone(T,1);
    
    return;
  }

  template<typename CoeffType>
  inline
  LindemannFalloff<CoeffType>::LindemannFalloff(const unsigned int nspec):n_spec(nspec)
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
