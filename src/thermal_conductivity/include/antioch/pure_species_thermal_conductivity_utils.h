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

#ifndef ANTIOCH_PURE_SPECIES_THERMAL_CONDUCTIVITY_UTILS_H
#define ANTIOCH_PURE_SPECIES_THERMAL_CONDUCTIVITY_UTILS_H

#include "antioch/antioch_asserts.h"
#include "antioch/pure_species_the_con_utils_decl.h"
#include "antioch/pure_species_thermal_conductivity_building.h"

namespace Antioch
{
   // getting tag
   template <typename ThermoEvaluator, typename CoeffType>
   struct physical_tag<PureSpeciesThermalConductivity<ThermoEvaluator,CoeffType> >:
        public physical_tag_base<PureSpeciesThermalConductivity<ThermoEvaluator,CoeffType> >
   {
      typedef pure_species_thermal_conductivity_tag type;
        // some models require specific initialization
        // (typically automatic initialization)
     typedef pure_species_thermal_conductivity_tag init_type;
     typedef pure_species_thermal_conductivity_tag set_type;
        // for operators, diffusion is special, see comment below
     typedef pure_species_thermal_conductivity_tag thermal_conductivity_type;
   };

   // physical set boolean
   template<typename ThermoEvaluator, typename CoeffType>
   struct is_physical_set<PureSpeciesThermalConductivity<ThermoEvaluator,CoeffType> >
   {
      static const bool value = true;
   };

   template <typename ThermoEvaluator, typename CoeffType>
   struct Initializer<PureSpeciesThermalConductivity<ThermoEvaluator,CoeffType>, pure_species_thermal_conductivity_tag >
  {
     Initializer(const ThermoEvaluator & th, const CoeffType  Z_298K_, const CoeffType  LJ_depth_):
                        t(th),Z_298K(Z_298K_),LJ_depth(LJ_depth_){}

     const ThermoEvaluator & t;
     const CoeffType       & Z_298K;
     const CoeffType       & LJ_depth;
  };

   // custom add
   template <typename Model, typename InitType>
   void physical_set_add(unsigned int s, typename SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, pure_species_thermal_conductivity_tag)
   {
       antioch_assert_less(s,set.size());
       antioch_assert(!set[s]);

       set[s] = new Model(init.t, init.Z_298K, init.LJ_depth);
   }

   // requires to void this
   template <typename Model, typename InitType>
   void physical_set_add(Model & set, const InitType & init, pure_species_thermal_conductivity_tag){}

   // custom reset
   template <typename Model, typename InitType>
   void physical_set_reset(unsigned int s, typename SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init,  pure_species_thermal_conductivity_tag)
   {
       antioch_assert_less(s,set.size());
       antioch_assert(set[s]);

       set[s]->reset_coeffs(init.t, init.Z_298K, init.M, init.LJ_depth);
   }

   // requires to void this
   template <typename Model, typename InitType>
   void physical_set_reset(Model & set, const InitType & init,  pure_species_thermal_conductivity_tag){}


   // operator
   template<typename Model, typename StateType>
   void physical_set_operator_thermal_conductivity(const Model & set, unsigned int s, const StateType & mu, const StateType & dss, const StateType & T, const StateType & rho, 
                                                   StateType & k, pure_species_thermal_conductivity_tag)
   {
        antioch_assert(!set.empty());

        k = (*set[s])(s,mu,T,rho,dss);
   }

   template<typename Model, typename StateType, typename VectorStateType>
   void physical_set_operator_thermal_conductivity(const Model & set, const VectorStateType & mu, const VectorStateType & dss, const StateType & T, const StateType & rho, 
                                                   VectorStateType & k, pure_species_thermal_conductivity_tag)
   {
      antioch_assert_equal_to(k.size(),set.size());

      for(unsigned int s = 0; s < k.size(); ++s)
      {
          k[s] = (*set[s])(s,mu[s],T,rho,dss[s]);
      }
   }

   template <typename ModelSet>
   void physical_set_initialize(ModelSet & mod, pure_species_thermal_conductivity_tag)
   {
      mod.set().resize(mod.mixture().n_species(),NULL);
      build_pure_species_thermal_conductivity(mod);
   }

}

#endif
