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

// Antioch
#include "antioch/metaprogramming_decl.h"
#include "antioch/antioch_asserts.h"

// C++
#include <map>

namespace Antioch
{

   template <typename Physics, class Enable = void>
   struct is_physical_set
   {
      static const bool value = false;
   };

   // more explicit name
   typedef numeric_library_tag default_physical_tag;
   struct default_physical_set_tag {};

// gives the tag of the model for any interesting
// action, base class so we don't rewrite everything
   template <typename Model, typename Enable = void>
   struct physical_tag_base
   {
        // tag
     typedef default_physical_tag type;
        // kind of set tag
     typedef typename if_else_type<is_physical_set<Model>::value,
                                     default_physical_set_tag,
                                     default_physical_tag
                                   >::type set_type;
        // some models require specific initialization
        // (typically automatic initialization)
     typedef typename if_else_type<is_physical_set<Model>::value,
                                     default_physical_set_tag,
                                     default_physical_tag
                                   >::type init_type;
        // some models require specific initialization
        // but not specific deletion
     typedef typename if_else_type<is_physical_set<Model>::value,
                                     default_physical_set_tag,
                                     default_physical_tag
                                   >::type del_type;
        // for operators, diffusion is special, see comment below
     typedef default_physical_tag viscosity_type;
     typedef default_physical_tag diffusion_species_type;
     typedef default_physical_tag diffusion_mixture_type;
     typedef default_physical_tag thermal_conductivity_type;
   };

   template <typename Model, typename Enable = void>
   struct physical_tag: public physical_tag_base<Model,Enable>
   {};

   /* Whether we have an equation
      for every species, or we need
      a set of species (typically a std::vector<Model*>).
      default is an equation
   */
   template <typename Model, bool B>
   struct SetOrEquation
   {
      typedef Model* type;
   };

   template <typename Model>
   struct SetOrEquation<Model,true>
   {
      typedef std::vector<Model*> type;
   };

  // might be required when std::vector<> can't hold all the
  // initialization
  template<typename Physics, typename Tag = default_physical_tag>
  struct Initializer{};

    // init
   template <typename ModelSet>
   void physical_set_initialize(ModelSet & mod, default_physical_tag );

   template <typename ModelSet>
   void physical_set_initialize(ModelSet & mod, default_physical_set_tag );

    // deleted
   template <typename Model>
   void physical_set_delete(Model & set, default_physical_tag);

   template <typename Model>
   void physical_set_delete(Model & set, default_physical_set_tag);

   // add general model
   template <typename Model, typename InitType>
   void physical_set_add(typename SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, default_physical_tag);

   template <typename Model, typename InitType>
   void physical_set_add(typename SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, default_physical_set_tag);

   // add species model
   template <typename Model, typename InitType>
   void physical_set_add(unsigned int s, Model & set, const InitType & init, default_physical_tag);

   template <typename Model, typename InitType>
   void physical_set_add(unsigned int s, Model & set, const InitType & init, default_physical_set_tag);

   // reset general model
   template <typename Model, typename InitType>
   void physical_set_reset(Model & set, const InitType & init, default_physical_tag);

   template <typename Model, typename InitType>
   void physical_set_reset(Model & set, const InitType & init, default_physical_set_tag);

   // reset species model
   template <typename Model, typename InitType>
   void physical_set_reset(unsigned int s, Model & set, const InitType & init, default_physical_tag);

   template <typename Model, typename InitType>
   void physical_set_reset(unsigned int s, Model & set, const InitType & init, default_physical_set_tag);

   // print
   template <typename Model>
   void physical_set_print(const Model & set, const std::map<unsigned int, std::string> & spec, std::ostream & out, default_physical_tag);

   template <typename Model>
   void physical_set_print(const Model & set, const std::map<unsigned int, std::string> & spec, std::ostream & out, default_physical_set_tag);

    // operators
    /*
        Each type of tag (viscosity, diffusion, thermal conduction)
        will enable two (and only two) of those methods.
        In the special case of diffusion, we need to disable
        the two methods that will have a non-default signature.
        This is because binary molecular model is indeed
        one more dimension above mixture and not constant
        Lewis. Thus they act differently, thus two tags
        as a strategy to disable.
     */

    ///// viscosity

        /// may or may not be used, though always
        /// defined
        /// should not be used when not concerned
   template<typename Model, typename StateType>
   void physical_set_operator_viscosity(const Model & set, unsigned int s, const StateType & T, StateType & k, default_physical_tag);

        /// this one may or may not be used,
        /// void by default to skip it in concerned cases
   template<typename Model, typename StateType, typename VectorStateType>
   void physical_set_operator_viscosity(const Model & set, const StateType & T, VectorStateType & mu, default_physical_tag);


    ///// thermal conduction

        /// may or may not be used, though always
        /// defined
        /// should not be used in not concerned
   template<typename Model, typename StateType>
   void physical_set_operator_thermal_conductivity(const Model & set, unsigned int s, const StateType & mu, 
                                                   const StateType & dss, const StateType & T, const StateType & rho,
                                                   StateType & k,
                                                    default_physical_tag);

        /// this one may or may not be used
        /// nothing is not needed
   template<typename Model, typename StateType, typename VectorStateType>
   void physical_set_operator_thermal_conductivity(const Model & set, const VectorStateType & mu, 
                                                   const VectorStateType & dss, const StateType & T, const StateType & rho, 
                                                    default_physical_tag);


    ///// diffusion

        /// this one may or may not be used,
        /// void by default to skip it in concerned cases
   template<typename Model, typename StateType, typename MatrixStateType>
   void physical_set_operator_diffusion(const Model & set, const StateType & T, const StateType & cTot,
                                     MatrixStateType & Ds, default_physical_tag);

        /// this one may or may not be used,
        /// void by default to skip it in concerned cases
   template<typename Model, typename StateType>
   void physical_set_operator_diffusion(unsigned int s, const Model & set, const StateType & T, const StateType & cTot,
                                         StateType & dss, default_physical_tag);

        /// this one may or may not be used,
        /// void by default to skip it in concerned cases
   template<typename Model, typename StateType>
   void physical_set_operator_diffusion(const Model & set, const StateType & rho, const StateType & cp, const StateType & k, StateType & ds, default_physical_tag);

        /// this one may or may not be used,
        /// void by default to skip it in concerned cases
   template<typename Model, typename StateType, typename VectorStateType>
   void physical_set_operator_diffusion(const Model & set, const StateType & rho, const VectorStateType & cp, const VectorStateType & k, VectorStateType & ds, default_physical_tag);

        /// To have Wilke rule inside Wilke object and not
        /// in diffusion object, avoid it by default to
        /// skip it in concerned cases
   template <typename Mixture, typename MatrixStateType, typename VectorStateType>
   void wilke_diffusion_rule(const Mixture & mixture, const VectorStateType & mass_fractions, const MatrixStateType & Ds, VectorStateType & ds, default_physical_tag);
}

#endif
