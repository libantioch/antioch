//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
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

#ifndef ANTIOCH_KINETICS_THEORY_THERMAL_CONDUCTIVITY_H
#define ANTIOCH_KINETICS_THEORY_THERMAL_CONDUCTIVITY_H

//Antioch
#include "antioch/antioch_assert.h"
#include "antioch/thermal_conductivity_enum.h"
#include "antioch/rotational_relaxation.h"
#include "antioch/metaprogramming_decl.h"

//C++


namespace Antioch{

  template <typename CoeffType>
  class KineticsTheoryThermalConductivity
  {
      public:

        KineticsTheoryThermalConductivity(const CoeffType & Z_298K, const CoeffType & M, const CoeffType & LJ_depth);
                                   

        ~KineticsTheoryThermalConductivity();


        template<typename StateType>
        const
        ANTIOCH_AUTO(StateType
        operator()(const StateType & T, const StateType &molar_concentration, //state
                                  const StateType & Cv_vib,const StateType & Cv_rot,const StateType & Cv_trans, //Cv
                                  const StateType & vis, const StateType & self_diffusion) const//other transports

        ANTIOCH_AUTOFUNC(StateType,this->thermal_conductivity(T, molar_concentration, Cv_vib, Cv_rot, Cv_trans, vis, self_diffusion))



        template<typename StateType>
        const StateType thermal_conductivity(const StateType & T, const StateType &molar_concentration, //state
                                             const StateType & Cv_vib,const StateType & Cv_rot,const StateType & Cv_trans, //Cv
                                             const StateType & vis, const StateType & self_diffusion) const; //other transports


        ThermalConductivityModel model() const {return _model;}

      private:

        /*! never ever use it*/
        KineticsTheoryThermalConductivity();

//small enough
        RotationalRelaxation<CoeffType> _rot;

        CoeffType                       _M;        // molar mass
        CoeffType                       _m;        // molecular mass
        CoeffType                       _LJ_depth; // Lennard-Jones depth

//constants
        const CoeffType five_over_two;
        const CoeffType five_over_three;
        const CoeffType two_over_pi;
        const CoeffType one;

        const ThermalConductivityModel _model;

  };

  template <typename CoeffType>
  inline
  KineticsTheoryThermalConductivity<CoeffType>::KineticsTheoryThermalConductivity(const CoeffType & Z_298K, const CoeffType & M, const CoeffType & LJ_depth):
        _rot(Z_298K),
        _M(M),
        _m(_M/Constants::Avogadro<StateType>()),
        _LJ_depth(LJ_depth),
        five_over_two(5.L/2.L),
        five_over_three(5.L/3.L),
        two_over_pi(2.L/Constants::pi<CoeffType>()),
        one(1.L),
        _model(PURE_SPECIES)
  {
      return;
  }

  template <typename CoeffType>
  inline
  KineticsTheoryThermalConductivity<CoeffType>::~KineticsTheoryThermalConductivity()
  {
     return;
  }

  template <typename CoeffType>
  template <typename StateType, typename VectorStateType >
  inline
  void KineticsTheoryThermalConductivity<CoeffType>::thermal_conductivity(const StateType & T, const StateType &molar_concentration, //state
                                                                   const StateType & Cv_vib, const StateType & Cv_rot, const StateType & Cv_trans, //Cv
                                                                   const StateType & vis, const StateType & self_diffusion) const //other transports
  {

// tmp data, compute once
    StateType rho_times_D_over_vis = molar_concentration * _M * self_diffusion / vis;
    StateType A                    = five_over_two - rho_times_D_over_vis ;
    StateType B                    = _rot(T,_LJ_depth) + 
                                     two_over_pi * (five_over_three * Cv_rot / Constants::R_universal<StateType> * StateType(1000.)  + // to SI
                                                    rho_times_D_over_vis ) ; 
//
    return  vis / _m * ( five_over_two   * (one - two_over_pi * (Cv_rot * A ) / (Cv_trans * B ) ) * Cv_trans // .. ( f_trans * Cv,trans + ..
                         + ( rho_times_D_over_vis  * (one + two_over_pi * A / B)) * Cv_rot // ... + f_rot * Cv_rot + .. 
                         + ( rho_times_D_over_vis) * Cv_vib //  ..f_vib * Cv_vib);
                       );   
  }

} //end namespace Antioch

#endif
