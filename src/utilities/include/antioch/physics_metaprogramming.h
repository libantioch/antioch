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
   template <typename ModelSet>
   void physical_set_initialize(ModelSet & mod, default_physical_tag ){mod.set() = NULL;}

   template <typename ModelSet>
   void physical_set_initialize(ModelSet & mod, default_physical_set_tag ){mod.set().resize(mod.mixture().n_species(),NULL);}

   // delete
   template <typename ModelSet>
   void physical_set_delete(ModelSet & mod, default_physical_tag) 
   {
      if(mod.set())delete mod.set();
      mod.set() = NULL;
   }

   template <typename ModelSet>
   void physical_set_delete(ModelSet & mod, default_physical_set_tag)
   {
      for(unsigned int s = 0; s < mod.mixture().n_species(); ++s)
      {
         if(mod.set()[s])delete mod.set()[s];
         mod.set()[s] = NULL;
      }
   }

   // add general model
   template <typename Model, typename InitType>
   void physical_set_add(typename SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, default_physical_tag)
   {
      antioch_assert(!set);
      set = new Model(init);
   }

   template <typename Model, typename InitType>
   void physical_set_add(typename SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, default_physical_set_tag){}

   // add species model
   template <typename Model, typename InitType>
   void physical_set_add(unsigned int s, typename SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, default_physical_tag){}

   template <typename Model, typename InitType>
   void physical_set_add(unsigned int s, typename SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, default_physical_set_tag)
   {
      antioch_assert_less(s,set.size());
      antioch_assert(!set[s]);

      set[s] = new Model(init);
   }

   // reset general model
   template <typename Model, typename InitType>
   void physical_set_reset(Model & set, const InitType & init, default_physical_tag)
   {
      set->reset_coeff(init);
   }

   template <typename Model, typename InitType>
   void physical_set_reset(Model & set, const InitType & init, default_physical_set_tag){}

   // reset species model
   template <typename Model, typename InitType>
   void physical_set_reset(unsigned int s, Model & set, const InitType & init, default_physical_tag){}

   template <typename Model, typename InitType>
   void physical_set_reset(unsigned int s, Model & set, const InitType & init, default_physical_set_tag)
   {
      set[s]->reset_coeffs(init);
   }

   // print
   template <typename Model>
   void physical_set_print(const Model & set, const std::map<unsigned int, std::string> & spec, std::ostream & out, default_physical_tag)
   {
      out << "Model used for species:" << std::endl;
      for(std::map<unsigned int, std::string>::const_iterator it = spec.begin(); it != spec.end(); ++it)
      {
          out << it->second << " "; 
      }
      out << std::endl << (*set) << std::endl;
   }

   template <typename Model>
   void physical_set_print(const Model & set, const std::map<unsigned int, std::string> & spec, std::ostream & out, default_physical_set_tag)
   {
      out << "Model per species: " << std::endl;

      for(unsigned int s = 0; s < spec.size(); ++s)
      {
          out << "\t" << spec.at(s) << ": " << *set[s] << std::endl;
      }
   }


/// operators

   /// viscosity
   // keep it this way of void it?
   template<typename Model, typename StateType>
   void physical_set_operator_viscosity(const Model & set, unsigned int s, const StateType & T, StateType & mu, default_physical_tag)
   {}

   // supposed to disappear and be optimized out if not needed
   template<typename Model, typename StateType, typename VectorStateType>
   void physical_set_operator_viscosity(const Model & set, const StateType & T, VectorStateType & mu, default_physical_tag)
   {}

   /// thermal conduction
   // keep it this way of void it?
   template<typename Model, typename StateType>
   void physical_set_operator_thermal_conductivity(const Model & set, unsigned int s, const StateType & mu, const StateType & dss, const StateType & T, const StateType & rho, 
                                                    StateType & k, default_physical_tag)
   {}

   // supposed to disappear and be optimized out if not needed
   template<typename Model, typename StateType, typename VectorStateType>
   void physical_set_operator_thermal_conductivity(const Model & set, const VectorStateType & mu, 
                                                   const VectorStateType & dss, const StateType & T, const StateType & rho, 
                                                    default_physical_tag)
   {}


   /// diffusion

   // supposed to disappear and be optimized out
   template<typename Model, typename StateType, typename MatrixStateType>
   void physical_set_operator_diffusion(const Model & set, const StateType & T, const StateType & cTot,
                                        MatrixStateType & Ds, default_physical_tag)
   {}

   // suppose to disappear and be optimized out
   template<typename Model, typename StateType>
   void physical_set_operator_diffusion(unsigned int s, const Model & set, const StateType & T, const StateType & cTot,
                                         StateType & dss, default_physical_tag)
   {}

   // suppose to disappear and be optimized out
   template<typename Model, typename StateType>
   void physical_set_operator_diffusion(const Model & set, const StateType & rho, const StateType & cp, const StateType & k, StateType & ds, default_physical_tag)
   {}

   // suppose to disappear and be optimized out
   template<typename Model, typename StateType, typename VectorStateType>
   void physical_set_operator_diffusion(const Model & set, const StateType & rho, const VectorStateType & cp, const VectorStateType & k, VectorStateType & ds, default_physical_tag)
   {}

   // suppose to disappear and be optimized out
   template <typename Mixture, typename MatrixStateType, typename VectorStateType>
   void wilke_diffusion_rule(const Mixture & mixture, const VectorStateType & mass_fractions, const MatrixStateType & Ds, VectorStateType & ds, default_physical_tag)
   {}

}

#endif
