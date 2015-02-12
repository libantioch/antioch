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

#ifndef ANTIOCH_THERMO_HANDLER_H
#define ANTIOCH_THERMO_HANDLER_H

#include "antioch/ideal_gas_micro_thermo.h"

namespace Antioch
{

  /* we want as macro values
        - cp, cp/R
        - cv, cv/R
        - h, h/(RT)
        - s, s/R
     as micro values
        - cv_trans, cp_trans/R
        - cv_rot, cp_rot/R
        - cv_vib, cp_vib/R

     Any other quantity will require to fetch
     the concerned thermo object. Those
     quantities are given per species, for
     mixture properties, fetch the concerned
     object.
  */
  template <typename CoeffType, typename MacroThermo, typename MicroThermo = IdealGasMicroThermo<MacroThermo,CoeffType> >
  class ThermoHandler
  {
      public:
        ThermoHandler(const MacroThermo & macro_thermo, const MicroThermo & micro_thermo);
        ~ThermoHandler();

// macro values

        //! Cp
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
           cp(const TempCache<StateType> & cache, unsigned int s) const
         ANTIOCH_AUTOFUNC(StateType, _macro_therm.cp(cache,s))

        //! Cp/R
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
           cp_over_R(const TempCache<StateType> & cache, unsigned int s) const
         ANTIOCH_AUTOFUNC(StateType, _macro_therm.cp_over_R(cache,s))

        //! Cv
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
           cv(const TempCache<StateType> & cache, unsigned int s) const
         ANTIOCH_AUTOFUNC(StateType, _macro_therm.cv(cache,s))

        //! Cv/R
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
           cv_over_R(const TempCache<StateType> & cache, unsigned int s) const
         ANTIOCH_AUTOFUNC(StateType, _macro_therm.cv_over_R(cache,s))

        //! h
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
           h(const TempCache<StateType> & cache, unsigned int s) const
         ANTIOCH_AUTOFUNC(StateType, _macro_therm.h(cache,s))

        //! h/(RT)
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
           h_over_RT(const TempCache<StateType> & cache, unsigned int s) const
         ANTIOCH_AUTOFUNC(StateType, _macro_therm.h_over_RT(cache,s))

        //! s
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
           s(const TempCache<StateType> & cache, unsigned int s) const
         ANTIOCH_AUTOFUNC(StateType, _macro_therm.h(cache,s))

        //! s/R
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
           s_over_R(const TempCache<StateType> & cache, unsigned int s) const
         ANTIOCH_AUTOFUNC(StateType, _macro_therm.s_over_R(cache,s))


// Micro thermo
        //! cv_vib
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
           cv_vib(const StateType & T, unsigned int s) const
         ANTIOCH_AUTOFUNC(StateType, _micro_therm.cv_vib(s,T))

        //! cv_vib/R
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
           cv_vib_over_R(const StateType & T, unsigned int s) const
         ANTIOCH_AUTOFUNC(StateType, _micro_therm.cv_vib_over_R(s,T))

        //! cv_rot
        const CoeffType cv_rot(unsigned int s)          const {return _micro_therm.cv_vib(s);}

        //! cv_rot/R
        const CoeffType cv_rot_over_R(unsigned int s)   const {return _micro_therm.cv_rot_over_R(s);}

        //! cv_trans
        const CoeffType cv_trans(unsigned int s)        const {return _micro_therm.cv_trans(s);}

        //! cv_trans/R
        const CoeffType cv_trans_over_R(unsigned int s) const {return _micro_therm.cv_trans_over_R(s);}


// objects

        const MacroThermo & macro_thermo() const;

        const MicroThermo & micro_thermo() const;


      private:
        // huhu, don't use it
        ThermoHandler();

        const MacroThermo & _macro_therm;
        const MicroThermo & _micro_therm;
  };


   template <typename CoeffType, typename MacroThermo, typename MicroThermo>
   inline
   ThermoHandler<CoeffType,MacroThermo, MicroThermo>::ThermoHandler(const MacroThermo & macro_thermo, const MicroThermo & micro_thermo):
        _macro_therm(macro_thermo),
        _micro_therm(micro_thermo)
   {
      return;
   }

   template <typename CoeffType, typename MacroThermo, typename MicroThermo>
   inline
   ThermoHandler<CoeffType, MacroThermo, MicroThermo>::~ThermoHandler()
   {
     return;
   }


   template <typename CoeffType, typename MacroThermo, typename MicroThermo>
   inline
   const MacroThermo & ThermoHandler<CoeffType, MacroThermo, MicroThermo>::macro_thermo() const
   {
      return _macro_therm;
   }
   

   template <typename CoeffType, typename MacroThermo, typename MicroThermo>
   inline
   const MicroThermo &  ThermoHandler<CoeffType, MacroThermo, MicroThermo>::micro_thermo() const
   {
      return _micro_therm;
   }

}

#endif
