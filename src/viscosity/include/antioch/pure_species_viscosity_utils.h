//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_PURE_SPECIES_VISCOSITY_UTILS_H
#define ANTIOCH_PURE_SPECIES_VISCOSITY_UTILS_H

#include "antioch/pure_species_viscosity_utils_decl.h"

namespace Antioch
{
   // getting tag
   template <typename CoeffType, typename Interpolator>
   struct physical_tag_type<PureSpeciesViscosity<CoeffType, Interpolator> >
   {
      typedef pure_species_viscosity_tag type;
   };

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
