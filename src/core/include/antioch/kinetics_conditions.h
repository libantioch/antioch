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

#ifndef ANTIOCH_KINETICS_CONDITIONS_H
#define ANTIOCH_KINETICS_CONDITIONS_H

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/particle_flux.h"
#include "antioch/temp_cache.h"

// C++
#include <vector>
#include <map>

namespace Antioch{

  /*! This class contains the conditions of the chemistry

      Conditions are:
        - Temperature(s)
        - particle flux
        - any other not yet implemented conditions (catalysis for instance)

      Pressure is given by the molecular composition of the mixture,
      we are in a gaz phase, so the state equation will give pressure (typically
      ideal gas).

      The idea is to avoid any copy of anything as much as we can, 
      we store pointers but deal only with references. Nothing
      belongs to this object, if there is any cleaning to do,
      it must be done elsewhere.

      Might be interesting to see about the "double" case
      for temperature: faster to copy instead of copying
      the reference, and it's expected to be the classic
      case.
   */
  template <typename StateType, 
            typename VectorStateType = std::vector<StateType> >
  class KineticsConditions
  {
        public:

          KineticsConditions(const StateType & temperature);
          ~KineticsConditions();

          void add_particle_flux(const ParticleFlux<VectorStateType> & pf, unsigned int nr);

          const StateType & T() const;

          //! returns the temperature T,
          //
          // \todo generalize
          const StateType & Tvib() const;

          const TempCache<StateType> & temp_cache() const;

          const ParticleFlux<VectorStateType> & particle_flux(int nr) const;

        private:

          KineticsConditions();

          TempCache<StateType> _temperature; 

        // pointer's not const, particle flux is
          std::map<unsigned int,ParticleFlux<VectorStateType> const * const > _map_pf; 

  };


  template <typename StateType, typename VectorStateType>
  inline
  KineticsConditions<StateType,VectorStateType>::KineticsConditions(const StateType & temperature):
        _temperature(temperature)
  {
    return;
  }

  template <typename StateType, typename VectorStateType>
  inline
  KineticsConditions<StateType,VectorStateType>::~KineticsConditions()
  {
    return;
  }

  template <typename StateType, typename VectorStateType>
  inline
  void KineticsConditions<StateType,VectorStateType>::add_particle_flux(const ParticleFlux<VectorStateType> & pf, unsigned int nr)
  {
     _map_pf.insert(std::make_pair(nr, &pf));
  }

  template <typename StateType, typename VectorStateType>
  inline
  const StateType & KineticsConditions<StateType,VectorStateType>::T() const
  {
     return _temperature.T;
  }

  template <typename StateType, typename VectorStateType>
  inline
  const StateType & KineticsConditions<StateType,VectorStateType>::Tvib() const
  {
     return _temperature.T;
  }

  template <typename StateType, typename VectorStateType>
  inline
  const ParticleFlux<VectorStateType> & KineticsConditions<StateType,VectorStateType>::particle_flux(int nr) const
  {
     antioch_assert(_map_pf.count(nr));
     return *(_map_pf.at(nr));
  }

  template <typename StateType, typename VectorStateType>
  inline
  const TempCache<StateType> & KineticsConditions<StateType,VectorStateType>::temp_cache() const
  {
     return _temperature;
  }

 // partial specialization for the return_type
  template <typename StateType, typename VectorStateType>
  struct return_type<KineticsConditions<StateType,VectorStateType> >
  {
     typedef StateType type;
  };

} //end namespace Antioch

#endif
