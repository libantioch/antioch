//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_PURE_SPECIES_VISCOSITY_UTILS_H
#define ANTIOCH_PURE_SPECIES_VISCOSITY_UTILS_H

#include "antioch/antioch_asserts.h"
#include "antioch/pure_species_viscosity_utils_decl.h"
#include "antioch/pure_species_viscosity_building.h"

namespace Antioch
{
   // getting tag
   template <typename CoeffType, typename Interpolator>
   struct physical_tag<PureSpeciesViscosity<CoeffType, Interpolator> >
   {
      typedef PureSpeciesViscosity<CoeffType, Interpolator> Model;
      typedef pure_species_viscosity_tag type;
        // kind of set tag
     typedef typename if_else_type<is_physical_set<Model>::value,
                                     default_physical_set_tag,
                                     default_physical_tag
                                   >::type set_type;
        // some models require specific initialization
        // (typically automatic initialization)
     typedef pure_species_viscosity_tag init_type;
        // some models require specific initialization
        // but not specific deletion
     typedef typename if_else_type<is_physical_set<Model>::value,
                                     default_physical_set_tag,
                                     default_physical_tag
                                   >::type del_type;
        // for operators, diffusion is special, see comment below
     typedef pure_species_viscosity_tag viscosity_type;
     typedef default_physical_tag       diffusion_species_type;
     typedef default_physical_tag       diffusion_mixture_type;
     typedef default_physical_tag       thermal_conductivity_type;
   };

   // physical set boolean
   template<typename CoeffType, typename Interpolator>
   struct is_physical_set<PureSpeciesViscosity<CoeffType, Interpolator> >
   {
      static const bool value = true;
   };

   // we can initialize without the user's help,
   // so we tag the physical_set_tag:
   // we need to define 
   // physical_set_initialize
   // physical_set_delete
   // physical_set_add
   // physical_set_reset
   // but not the default ones that are not concerned
   template <typename CoeffType, typename Interpolator>
   struct physical_set_tag<PureSpeciesViscosity<CoeffType,Interpolator> >
   {
     typedef typename pure_species_viscosity_tag type;
   };

   // we can initialize without the user's help
   template <typename ModelSet>
   void physical_set_initialize(ModelSet & mod, pure_species_viscosity_tag )
   {
      mod.set().resize(mod.mixture().n_species(),NULL);
      build_pure_species_viscosity(mod);
   }

   template<typename Model, typename StateType>
   void physical_set_operator_viscosity(const Model & set, unsigned int s, const StateType & T, StateType & mu, pure_species_viscosity_tag)
   {
      mu = (*set[s])(T);
   }

   template<typename Model, typename StateType, typename VectorStateType>
   void physical_set_operator_viscosity(const Model & set, const StateType & T, VectorStateType & mu, sutherland_viscosity_tag)
   {
      antioch_assert_equal_to(mu.size(), set.size());

      for(unsigned int s = 0; s < mu.size(); s++)
      {
          mu[s] = (*set[s])(T);
      }
   }

}


#endif
