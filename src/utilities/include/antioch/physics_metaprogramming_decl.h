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

#ifndef ANTIOCH_PHYSICS_METAPROGRAMMING_DECL_H
#define ANTIOCH_PHYSICS_METAPROGRAMMING_DECL_H

#include "antioch/metaprogramming_decl.h"
#include "antioch/antioch_asserts.h"

namespace Antioch
{

   template <typename Physics, class Enable = void>
   struct is_physical_set
   {
      const bool value = false;
   };

   struct default_physical_tag {};
   struct physical_set_type_tag {};

// gives the tag of the model
   template <typename Model, typename Enable = void>
   struct physical_tag
   {
     typedef default_physical_tag type;
     typedef default_physical_tag second_type;
   };

// tag for set of models or just an equation
   template <typename Model, typename Enable = void>
   struct physical_set_tag
   {
     typedef typename if_else_type<is_physical_set<Model>,
                                   default_physical_tag, 
                                    physical_set_tag>::type type;
   };

   /* Wheter we have an equation
      for every species, or we need
      a set of species (typically a std::vector<Model*>).
      default is an equation (as of now, the most common)
   */
   template <typename Model, bool B>
   struct SetOrEquation
   {
      typedef Model * type;
   };

   template <typename Model>
   struct SetOrEquation<Model,true>
   {
      typedef std::vector<Model> type;
   };

  // might be required when std::vector<> can't hold all the
  // initialization
  template<typename Physics>
  struct Initializer;

    // init
   template <typename Model>
   void physical_set_initialize(Model & mod, default_physical_tag );

   template <typename Model>
   void physical_set_initialize(Model & mod, physical_set_type_tag );

    // deleted
   template <typename Model>
   void physical_set_delete(Model & set, default_physical_tag);

   template <typename Model>
   void physical_set_delete(Model & set, physical_set_tag);

   // add general model
   template <typename Model, typename InitType>
   void physical_set_add(SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, default_physical_tag);

   template <typename Model, typename InitType>
   void physical_set_add(SetOrEquation<Model,is_physical_set<Model>::value>::type  & set, const InitType & init, physical_set_tag);

   // add species model
   template <typename Model, typename InitType>
   void physical_set_add(unsigned int s, SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, default_physical_tag);

   template <typename Model, typename InitType>
   void physical_set_add(unsigned int s, SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, physical_set_tag);

   // reset general model
   template <typename Model, typename InitType>
   void physical_set_reset(Model & set, const InitType & init, default_physical_tag);

   template <typename Model, typename InitType>
   void physical_set_reset(Model & set, const InitType & init, physical_set_tag);

   // reset species model
   template <typename Model, typename InitType>
   void physical_set_reset(unsigned int s, Model & set, const InitType & init, default_physical_tag);

   template <typename Model, typename InitType>
   void physical_set_reset(unsigned int s, Model & set, const InitType & init, physical_set_tag);

    // operators
        /// this one returns a value, should always be
        /// ok to use, so error if not used correctly
   template<typename Model, typename StateType>
   ANTIOCH_AUTO(StateType) physical_set_operator_viscosity(const Model & set, unsigned int s, const StateType & T, default_physical_tag);

        /// this one may or may not be used
        /// nothing is not needed
   template<typename Model, typename StateType, typename VectorStateType>
   void physical_set_operator_thermal_conductivity(const Model & set, unsigned int s, const StateType & mu, const StateType & dss, const StateType & T, const StateType & rho, 
                                                   StateType & k, default_physical_tag);
        /// this one may or may not be used,
        /// void by default to skip it in concerned cases
   template<typename Model, typename StateType, typename MatrixStateType>
   void physical_set_operator_diffusion_comes_first(const Model & set, const StateType & T, const StateType & cTot,
                                     MatrixStateType & Ds, default_physical_tag);

        /// this one may or may not be used,
        /// void by default to skip it in concerned cases
   template<typename Model, typename StateType>
   void physical_set_operator_diffusion_comes_last(const Model & set, const StateType & cp, const StateType & k, StateType & ds, default_physical_tag);

        /// To have Wilke rule inside Wilke object and not
        /// in diffusion object, avoid it by default to
        /// skip it in concerned cases
   template <typename MatrixStateType, typename VectorStateType>
   void wilke_diffusion_rule(const MatrixStateType & Ds, VectorStateType & ds,default_physical_tag);
}

#endif
