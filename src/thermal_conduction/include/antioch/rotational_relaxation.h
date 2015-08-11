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

#ifndef ANTIOCH_ROTATIONAL_RELAXATION_H
#define ANTIOCH_ROTATIONAL_RELAXATION_H

// Antioch
#include "antioch/metaprogramming_decl.h"
#include "antioch/math_constants.h"
#include "antioch/cmath_shims.h"

// C++

namespace Antioch
{
  template <typename CoeffType>
  class RotationalRelaxation
  {
      public:
        RotationalRelaxation( CoeffType z_298, CoeffType eps_kb);
        ~RotationalRelaxation(){};

        void reset_coeffs( CoeffType rot, CoeffType depth);

        // \todo
        // cache _eps_kb / 298
        template <typename StateType>
        ANTIOCH_AUTO(StateType)
          operator()(const StateType & T) const
        ANTIOCH_AUTOFUNC(StateType, _z_298 * this->F(_eps_kb / 298) / this->F(StateType(_eps_kb/T)))

        //!
        CoeffType Z_298() const
        { return _z_298; }

        //!
        CoeffType eps_over_kb() const
        { return _eps_kb; }

      private:

        RotationalRelaxation();

        template <typename StateType>
        ANTIOCH_AUTO(StateType)
          F(const StateType & eps_T) const
        ANTIOCH_AUTOFUNC(StateType, _one + _pi32_2 * ant_sqrt(eps_T) + _pi2_4_plus_2 * eps_T + _pi32 * ant_pow(eps_T,CoeffType(1.5))) //unambiguate pow() for vexcl

        CoeffType _z_298;
        CoeffType _eps_kb;
        const CoeffType _one;
        const CoeffType _pi32_2;
        const CoeffType _pi2_4_plus_2;
        const CoeffType _pi32;
  };

  template <typename CoeffType>
  RotationalRelaxation<CoeffType>::RotationalRelaxation( CoeffType z_298, CoeffType eps_kb):
                _z_298(z_298),
                _eps_kb(eps_kb),
                _one(1),
                _pi32_2(ant_pow(Constants::pi<CoeffType>(),CoeffType(1.5)) / 2),
                _pi2_4_plus_2(Constants::pi<CoeffType>() * Constants::pi<CoeffType>() / 4 + 2),
                _pi32(ant_pow(Constants::pi<CoeffType>(),CoeffType(1.5)))
  {
     return;
  }

  template <typename CoeffType>
  inline
  void RotationalRelaxation<CoeffType>::reset_coeffs( CoeffType rot, CoeffType depth)
  {
     _z_298 = rot;
     _eps_kb = depth;
  }

}

#endif
