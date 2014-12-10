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

#ifndef ANTIOCH_PHYSICAL_SET_H
#define ANTIOCH_PHYSICAL_SET_H

// Antioch
#include "antioch/physics_metaprogramming_decl.h"

// C++
#include <iostream>

namespace Antioch
{
  template<typename Physics, typename Mixture>
  class  PhysicalSet
  {
      public:
// always usable wherever
        typedef Physics model;

        PhysicalSet(const Mixture & mix);

        ~PhysicalSet();

        const Mixture & mixture() const;

        // writable reference to the set
        typename SetOrEquation<Physics, is_physical_set<Physics>::value>::type & set();

        // const reference to the set
        const typename SetOrEquation<Physics, is_physical_set<Physics>::value>::type & set() const;

        template <typename InitType>
        void add_model(const std::string & species_name, const InitType & initMe);

        template <typename InitType>
        void add_model(const InitType & initMe);

        template <typename InitType>
        void reset_model(unsigned int s, const InitType & initMe);

        template <typename InitType>
        void reset_model(const InitType & initMe);

        // viscosity one species
        template<typename StateType>
        void operator()(unsigned int s, const StateType & T, StateType & mu) const;

        // viscosity full set
        template<typename StateType, typename VectorStateType>
        void operator()(const StateType & T, VectorStateType & mu) const;

        // diffusion full set (set level, matrix)
        template<typename StateType, typename MatrixStateType>
        void operator()(const StateType & T, const StateType & cTot, MatrixStateType & Ds) const;

        // diffusion self coefficient
        template<typename StateType>
        void operator()(unsigned int s, const StateType & T, const StateType & cTot, StateType & dss) const;

        // diffusion species (mixture level, scalar)
        template<typename StateType>
        void operator()(const StateType & rho, const StateType & cp, const StateType & k, StateType & ds) const;

        // diffusion full set (mixture level, species)
        template<typename StateType, typename VectorStateType>
        void operator()(const StateType & rho, const VectorStateType & cp, const VectorStateType & k, VectorStateType & ds) const;

        // thermal conduction one species
        template<typename StateType>
        void operator()(unsigned int s, const StateType & mu, const StateType & dss, const StateType & T, const StateType & rho, StateType & k) const;

        // thermal conduction full set
        template<typename StateType, typename VectorStateType>
        void operator()(const VectorStateType & mu, const VectorStateType & dss, const StateType & T, const StateType & rho, VectorStateType & k) const;

        void print(std::ostream & out = std::cout) const;

        friend std::ostream & operator<< (std::ostream & out, const PhysicalSet<Physics,Mixture> & physical_set )
        {
            physical_set.print(out);
            return out;
        }
        

      private:
        const Mixture & _mixture;
        typename SetOrEquation<Physics,is_physical_set<Physics>::value>::type _set;
  };

  template<typename Physics, typename Mixture>
  inline
  PhysicalSet<Physics,Mixture>::PhysicalSet(const Mixture & mix):_mixture(mix)
  {
    physical_set_initialize(*this, typename physical_tag<Physics>::init_type ());
  }

  template<typename Physics, typename Mixture>
  inline
  PhysicalSet<Physics,Mixture>::~PhysicalSet()
  {
    physical_set_delete(*this, typename physical_tag<Physics>::del_type());
  }

  template<typename Physics, typename Mixture>
  inline
  const Mixture & PhysicalSet<Physics,Mixture>::mixture() const
  {
    return _mixture;
  }

        // writable reference to the set
  template<typename Physics, typename Mixture>
  inline
  typename SetOrEquation<Physics, is_physical_set<Physics>::value>::type & PhysicalSet<Physics,Mixture>::set()
  {
    return _set;
  }

        // const reference to the set
  template<typename Physics, typename Mixture>
  inline
  const typename SetOrEquation<Physics, is_physical_set<Physics>::value>::type & PhysicalSet<Physics,Mixture>::set() const 
  {
    return _set;
  }

  template<typename Physics, typename Mixture>
  template <typename InitType>
  inline
  void PhysicalSet<Physics,Mixture>::add_model(const std::string & species_name, const InitType & initMe)
  {
           
     antioch_assert( _mixture.species_name_map().find(species_name) !=
                     _mixture.species_name_map().end() );

     unsigned int s = _mixture.species_name_map().find(species_name)->second;

     physical_set_add<Physics>(s, _set, initMe, typename physical_tag<Physics>::set_type());
  }

  template<typename Physics, typename Mixture>
  template <typename InitType>
  inline
  void PhysicalSet<Physics,Mixture>::add_model(const InitType & initMe)
  {
    physical_set_add<Physics>(_set, initMe, typename physical_tag<Physics>::set_type());
  }

  template<typename Physics, typename Mixture>
  template <typename InitType>
  inline
  void PhysicalSet<Physics,Mixture>::reset_model(unsigned int s, const InitType & initMe)
  {
    physical_set_reset(s, _set, initMe, typename physical_tag<Physics>::set_type());
  }

  template<typename Physics, typename Mixture>
  template <typename InitType>
  inline
  void PhysicalSet<Physics,Mixture>::reset_model(const InitType & initMe)
  {
    physical_set_reset(_set, initMe, typename physical_tag<Physics>::set_type());
  }

  template<typename Physics, typename Mixture>
  inline
  void PhysicalSet<Physics,Mixture>::print(std::ostream & out) const
  {
    physical_set_print(_set, _mixture.species_inverse_name_map(), out, typename physical_tag<Physics>::set_type());  
  }

  // viscosity one
  template<typename Physics, typename Mixture>
  template<typename StateType>
  inline
  void PhysicalSet<Physics,Mixture>::operator()(unsigned int s, const StateType & T, StateType & mu) const
  {
     physical_set_operator_viscosity(_set,s,T, mu, typename physical_tag<Physics>::viscosity_type());
  }

  // viscosity full
  template<typename Physics, typename Mixture>
  template<typename StateType, typename VectorStateType >
  inline
  void PhysicalSet<Physics,Mixture>::operator()(const StateType & T, VectorStateType &mu) const
  {
      physical_set_operator_viscosity(_set,T,mu, typename physical_tag<Physics>::viscosity_type());
  }

  // diffusion full
  template<typename Physics, typename Mixture>
  template<typename StateType, typename MatrixStateType>
  inline
  void PhysicalSet<Physics,Mixture>::operator()(const StateType & T, const StateType & cTot, MatrixStateType & Ds) const 
  {
    physical_set_operator_diffusion(_set, T, cTot, Ds, typename physical_tag<Physics>::diffusion_species_type());
  }

  // diffusion self-diffusion
  template<typename Physics, typename Mixture>
  template<typename StateType>
  inline
  void PhysicalSet<Physics,Mixture>::operator()(unsigned int s, const StateType & T, const StateType & cTot, StateType & dss) const 
  {
    physical_set_operator_diffusion(s,_set, T, cTot, dss, typename physical_tag<Physics>::diffusion_species_type());
  }

  // diffusion one, mixture level
  template<typename Physics, typename Mixture>
  template<typename StateType>
  inline
  void PhysicalSet<Physics,Mixture>::operator()(const StateType & rho, const StateType & cp, const StateType & k, StateType & ds) const
  {
    physical_set_operator_diffusion(_set, rho, cp, k, ds, typename physical_tag<Physics>::diffusion_mixture_type());
  }

  // diffusion full, mixture level
  template<typename Physics, typename Mixture>
  template<typename StateType, typename VectorStateType>
  inline
  void PhysicalSet<Physics,Mixture>::operator()(const StateType & rho, const VectorStateType & cp, const VectorStateType & k, VectorStateType & ds) const
  {
    physical_set_operator_diffusion(_set, rho, cp, k, ds, typename physical_tag<Physics>::diffusion_mixture_type());
  }

  // thermal conductivity one
  template<typename Physics, typename Mixture>
  template<typename StateType>
  inline
  void PhysicalSet<Physics,Mixture>::operator()(unsigned int s, const StateType & mu, const StateType & dss, const StateType & T, const StateType & rho, StateType & k) const
  {
    physical_set_operator_thermal_conductivity(_set, s, mu, dss, T, rho, k, typename physical_tag<Physics>::thermal_conductivity_type());
  }

  // thermal conductivity full
  template<typename Physics, typename Mixture>
  template<typename StateType, typename VectorStateType>
  inline
  void PhysicalSet<Physics,Mixture>::operator()(const VectorStateType & mu, const VectorStateType & dss, const StateType & T, const StateType & rho, VectorStateType & k) const
  {
    physical_set_operator_thermal_conductivity(_set, mu, dss, T, rho, k, typename physical_tag<Physics>::thermal_conductivity_type());
  }
}

#endif
