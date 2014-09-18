//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_PURE_SPECIES_VISCOSITY_UTILS_H
#define ANTIOCH_PURE_SPECIES_VISCOSITY_UTILS_H

#include "antioch/pure_species_viscosity_utils_decl.h"
#include "antioch/pure_species_viscosity_building.h"

namespace Antioch
{
   // getting tag
   template <typename CoeffType, typename Interpolator>
   struct physical_tag_type<PureSpeciesViscosity<CoeffType, Interpolator> >
   {
      typedef pure_species_viscosity_tag type;
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
      build_pure_species_viscosity(mod);
   }

   template <typename ModelSet>
   void physical_set_delete(ModelSet & mod, pure_species_viscosity_tag )
   {}

   template <typename Model, typename InitType>
   void physical_set_add(unsigned int s, SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, pure_species_viscosity_tag)
   {
     antioch_assert(!set[s]);
     set[s] = new Model(init);
   }

   template <typename Model, typename InitType>
   void physical_set_rest(unsigned int s, SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, pure_species_viscosity_tag)
   {
     set[s]->reset_coeffs(init);
   }

   // physical set boolean
   template<typename CoeffType, typename Interpolator>
   struct is_physical_set<PureSpeciesViscosity<CoeffType, Interpolator> >
   {
      static const bool value = true;
   };

// Pure Species
   template<typename Model, typename StateType, typename VectorStateType>
   ANTIOCH_AUTO(StateType) 
        physical_set_first_operator(const Model & set, unsigned int s, const StateType & T, pure_species_viscosity_tag)
   ANTIOCH_AUTOFUNC(StateType,set[s](T))

}


#endif
