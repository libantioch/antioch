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

#ifndef ANTIOCH_THERMO_EVALUATOR_H
#define ANTIOCH_THERMO_EVALUATOR_H

#include "antioch/ideal_gas_internal_thermo.h"

namespace Antioch
{

  /* we want as external values
        - cp, cp/R
        - cv, cv/R
        - h, h/(RT)
        - s, s/R
     as internal values
        - cv_trans, cp_trans/R
        - cv_rot, cp_rot/R
        - cv_vib, cp_vib/R

     Any other quantity will require to fetch
     the concerned thermo object. Those
     quantities are given per species, for
     mixture properties, fetch the concerned
     object.
  */
  template <typename CoeffType, typename ExternalThermo, typename InternalThermo = IdealGasInternalThermo<ExternalThermo,CoeffType> >
  class ThermoEvaluator
  {
      public:
        ThermoEvaluator(const ExternalThermo & external_thermo, const InternalThermo & internal_thermo);
        ~ThermoEvaluator();

// external values

        //! Cp
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
           cp(const TempCache<StateType> & cache, unsigned int s) const
         ANTIOCH_AUTOFUNC(StateType, _ext_therm.cp(cache,s))

        //! Cp/R
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
           cp_over_R(const TempCache<StateType> & cache, unsigned int s) const
         ANTIOCH_AUTOFUNC(StateType, _ext_therm.cp_over_R(cache,s))

        //! Cv
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
           cv(const TempCache<StateType> & cache, unsigned int s) const
         ANTIOCH_AUTOFUNC(StateType, _ext_therm.cv(cache,s))

        //! Cv/R
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
           cv_over_R(const TempCache<StateType> & cache, unsigned int s) const
         ANTIOCH_AUTOFUNC(StateType, _ext_therm.cv_over_R(cache,s))

        //! h
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
           h(const TempCache<StateType> & cache, unsigned int s) const
         ANTIOCH_AUTOFUNC(StateType, _ext_therm.h(cache,s))

        //! h/(RT)
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
           h_over_RT(const TempCache<StateType> & cache, unsigned int s) const
         ANTIOCH_AUTOFUNC(StateType, _ext_therm.h_over_RT(cache,s))

        //! s
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
           s(const TempCache<StateType> & cache, unsigned int s) const
         ANTIOCH_AUTOFUNC(StateType, _ext_therm.h(cache,s))

        //! s/R
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
           s_over_R(const TempCache<StateType> & cache, unsigned int s) const
         ANTIOCH_AUTOFUNC(StateType, _ext_therm.s_over_R(cache,s))


// Internal thermo
        //! cv_vib
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
           cv_vib(const StateType & T, unsigned int s) const
         ANTIOCH_AUTOFUNC(StateType, _int_therm.cv_vib(s,T))

        //! cv_vib/R
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
           cv_vib_over_R(const StateType & T, unsigned int s) const
         ANTIOCH_AUTOFUNC(StateType, _int_therm.cv_vib_over_R(s,T))

        //! cv_rot
        template <typename StateType>
        const CoeffType cv_rot(unsigned int s) const {return _int_therm.cv_vib(s);}

        //! cv_rot/R
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
           cv_rot_over_R(unsigned int s) const 
        ANTIOCH_AUTOFUNC(_int_therm.cv_rot_over_R(s,T))

        //! cv_trans
        template <typename StateType>
        const CoeffType cv_trans(unsigned int s) const {return _int_therm.cv_trans(s);}

        //! cv_trans/R
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
          cv_trans_over_R(unsigned int s) const 
        ANTIOCH_AUTOFUNC(_int_therm.cv_trans_over_R(s,T))


// objects

        const ExternalThermo & external_thermo() const;

        const InternalThermo & internal_thermo() const;


      private:
        // huhu, don't use it
        ThermoEvaluator();

        const ExternalThermo & _ext_therm;
        const InternalThermo & _int_therm;
  };


   template <typename ExternalThermo, typename InternalThermo>
   inline
   ThermoEvaluator<ExternalThermo, InternalThermo>::ThermoEvaluator(const ExternalThermo & external_thermo, const InternalThermo & internal_thermo):
        _ext_therm(external_thermo),
        _int_therm(internal_thermo)
   {
      return;
   }

   template <typename ExternalThermo, typename InternalThermo>
   inline
   ThermoEvaluator<ExternalThermo, InternalThermo>::~ThermoEvaluator()
   {
     return;
   }


   template <typename ExternalThermo, typename InternalThermo>
   inline
   const ExternalThermo & ThermoEvaluator<ExternalThermo, InternalThermo>::external_thermo() const
   {
      return _ext_thermo;
   }
   

   const InternalThermo &  ThermoEvaluator<ExternalThermo, InternalThermo>::internal_thermo() const
   {
      return _int_thermo;
   }

}

#endif
