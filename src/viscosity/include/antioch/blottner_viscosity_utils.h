//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_BLOTTNER_VISCOSITY_UTILS_H
#define ANTIOCH_BLOTTNER_VISCOSITY_UTILS_H

#include "antioch/blottner_viscosity_utils_decl.h"

namespace Antioch
{
// Blottner
   template <typename CoeffType>
   struct physical_tag_type<BlottnerViscosity<CoeffType> >
   {
      typedef blottner_viscosity_tag type;
   };

   // physical set boolean
   template<typename CoeffType>
   struct is_physical_set<BlottnerViscosity<CoeffType> >
   {
      static const bool value = true;
   };

   template<typename Model, typename StateType, typename VectorStateType>
   ANTIOCH_AUTO(StateType) 
        physical_set_operator_viscosity(const Model & set, unsigned int s, const StateType & T, blottner_viscosity_tag)
   ANTIOCH_AUTOFUNC(StateType,set[s](T))
}


#endif
