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

#ifndef ANTIOCH_PHYSICS_METAPROGRAMMING_H
#define ANTIOCH_PHYSICS_METAPROGRAMMING_H

#include "antioch/physics_metaprogramming_decl.h"
#include "antioch/metaprogramming_decl.h"
#include "antioch/antioch_asserts.h"

namespace Antioch
{
   // init
   template <typename Model>
   void physical_set_initialize(Model & mod, default_physical_tag ){mod = NULL;}

   template <typename Model>
   void physical_set_initialize(Model & mod, physical_set_type_tag ){}

   // delete
   template <typename Model>
   void physical_set_delete(Model & set, default_physical_tag) {delete set;}

   template <typename Model>
   void physical_set_delete(Model & set, physical_set_tag) {}

   // add general model
   template <typename Model, typename InitType>
   void physical_set_add(SetOrEquation<Model,is_physical_set<Model>::value>::type  & set, const InitType & init, default_physical_tag)
   {
      set = new Model(init);
   }

   template <typename Model, typename InitType>
   void physical_set_add(SetOrEquation<Model,is_physical_set<Model>::value>::type  & set, const InitType & init, physical_set_tag){}

   // add species model
   template <typename Model, typename InitType>
   void physical_set_add(unsigned int s, SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, default_physical_tag){}

   template <typename Model, typename InitType>
   void physical_set_add(unsigned int s, SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, physical_set_tag)
   {
      set.push_back(Model(init));
   }

   // reset general model
   template <typename Model, typename InitType>
   void physical_set_reset(Model & set, const InitType & init, default_physical_tag)
   {
      set->reset_coeff(init);
   }

   template <typename Model, typename InitType>
   void physical_set_reset(Model & set, const InitType & init, physical_set_tag){}

   // reset species model
   template <typename Model, typename InitType>
   void physical_set_reset(unsigned int s, Model & set, const InitType & init, default_physical_tag){}

   template <typename Model, typename InitType>
   void physical_set_reset(unsigned int s, Model & set, const InitType & init, physical_set_tag)
   {
      set[s].reset_coeffs(init);
   }

   template<typename Model, typename StateType>
   ANTIOCH_AUTO(StateType)
        physical_set_operator_viscosity(const Model & set, unsigned int s, const StateType & T, default_physical_tag)
   {
        antioch_error();// something's undefined
        return -1;
   }

   // suppose to disappear and be optimized out
   template<typename Model, typename StateType, typename VectorStateType>
   void physical_set_operator_thermal_conductivity(const Model & set, unsigned int s, const StateType & mu, const StateType & dss, const StateType & T, const StateType & rho, 
                                                   StateType & k, default_physical_tag)
   {}

   // suppose to disappear and be optimized out
   template<typename Model, typename StateType, typename MatrixStateType>
   void physical_set_operator_diffusion_comes_first(const Model & set, const StateType & T, const StateType & cTot,
                                     MatrixStateType & Ds, default_physical_tag)
   {}

   // suppose to disappear and be optimized out
   template<typename Model, typename StateType>
   void physical_set_operator_diffusion_comes_last(const Model & set, const StateType & cp, const StateType & k, StateType & ds, default_physical_tag)
   {}

   // suppose to disappear and be optimized out
   template <typename MatrixStateType, typename VectorStateType>
   void wilke_diffusion_rule(const MatrixStateType & Ds, VectorStateType & ds,default_physical_tag)
   {}

}

#endif
