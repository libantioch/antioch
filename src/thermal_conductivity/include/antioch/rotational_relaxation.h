//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_ROTATIONAL_RELAXATION_H
#define ANTIOCH_ROTATIONAL_RELAXATION_H

// Antioch
#include "antioch/math_constants.h"

// C++

namespace Antioch
{
  template <typename CoeffType>
  class RotationalRelaxation
  {
      public:
        RotationalRelaxation():
                _one(1.),
                _pi32_2(ant_pow(Constants::pi<CoeffType>(),1.5)),
                _pi2_4_plus_2(Constants::pi<CoeffType>() * Constants::pi<CoeffType>() / 4.L + 2.L),
                _pi32(ant_pow(Constants::pi<CoeffType>(),1.5))
                {}
        ~RotationalRelaxation(){}

        template <typename StateType>
        ANTIOCH_AUTO(StateType)
          operator()(const StateType & eps_T)
        ANTIOCH_AUTOFUNC(StateType, _one + _pi32_2 * ant_sqrt(eps_T) + _pi2_4_plus_2 * eps_T + _pi32 * ant_pow(eps_T,1.5))

      private:
        const CoeffType _one;
        const CoeffType _pi32_2;
        const CoeffType _pi2_4_plus_2;
        const CoeffType _pi32
  };

}

#endif
