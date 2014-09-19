//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_BLOTTNER_VISCOSITY_UTILS_H
#define ANTIOCH_BLOTTNER_VISCOSITY_UTILS_H

#include "antioch/antioch_asserts.h"
#include "antioch/blottner_viscosity_utils_decl.h"

namespace Antioch
{
// Blottner
   template <typename CoeffType>
   struct physical_tag<BlottnerViscosity<CoeffType> >:
        public physical_tag_base<BlottnerViscosity<CoeffType> >
   {
      typedef blottner_viscosity_tag type;
        // for operators, diffusion is special, see comment below
     typedef blottner_viscosity_tag viscosity_type;
   };

   // physical set boolean
   template<typename CoeffType>
   struct is_physical_set<BlottnerViscosity<CoeffType> >
   {
      static const bool value = true;
   };

   template<typename Model, typename StateType>
   void physical_set_operator_viscosity(const Model & set, unsigned int s, const StateType & T, StateType & mu, blottner_viscosity_tag)
   {
       mu = (*set[s])(T);
   }

   template<typename Model, typename StateType, typename VectorStateType>
   void physical_set_operator_viscosity(const Model & set, const StateType & T, VectorStateType & mu, blottner_viscosity_tag)
   {
      antioch_assert_equal_to(mu.size(), set.size());

      for(unsigned int s = 0; s < mu.size(); s++)
      {
          mu[s] = (*set[s])(T);
      }
   }
}


#endif
