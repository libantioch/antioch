//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_SUTHERLAND_VISCOSITY_UTILS_H
#define ANTIOCH_SUTHERLAND_VISCOSITY_UTILS_H

#include "antioch/sutherland_viscosity_utils_decl.h"

namespace Antioch
{
// Sutherland
   // getting tag
   template <typename CoeffType>
   struct physical_tag_type<SutherlandViscosity<CoeffType> >
   {
      typedef sutherland_viscosity_tag type;
   };

   // physical set boolean
   template<typename CoeffType>
   struct is_physical_set<SutherlandViscosity<CoeffType> >
   {
      static const bool value = true;
   };

   template<typename Model, typename StateType, typename VectorStateType>
   ANTIOCH_AUTO(StateType) 
        physical_set_operator_viscosity(const Model & set, unsigned int s, const StateType & T, sutherland_viscosity_tag)
   ANTIOCH_AUTOFUNC(StateType,set[s](T))
}


#endif
