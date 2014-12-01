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

#ifndef ANTIOCH_IDEAL_GAS_INTERNAL_THERMO_H
#define ANTIOCH_IDEAL_GAS_INTERNAL_THERMO_H

namespace Antioch
{

  template <typename ExternalThermo, typename CoeffType = double>
  class IdealGasInternalThermo
  {
       public:
         IdealGasInternalThermo(const ExternalThermo & ext_thermo, const ChemicalMixture<CoeffType> & chem_mix);
         ~IdealGasInternalThermo();

        //! cv_vib
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
           cv_vib(const StateType & T, unsigned int s) const
         ANTIOCH_AUTOFUNC(StateType, _ext_therm.cv(s,T) - this->cv_tr(s))

        //! cv_vib/R
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
           cv_vib_over_R(const StateType & T, unsigned int s) const
         ANTIOCH_AUTOFUNC(StateType, _ext_therm.cv_vib_over_R(s,T) - this->cv_tr_over_R(s))

        //! cv_rot
        const CoeffType cv_rot(unsigned int s) const;

        //! cv_rot/R
        template <typename StateType>
        const CoeffType cv_rot_over_R(unsigned int s) const;

        //! cv_trans
        template <typename StateType>
        const CoeffType cv_trans(unsigned int s) const;

        //! cv_trans/R
        template <typename StateType>
        const CoeffType cv_trans_over_R(unsigned int s) const;

        //! cv_trans-rot
        template <typename StateType>
        const CoeffType cv_tr(unsigned int s) const;

        //! cv_trans_rot/R
        template <typename StateType>
        const CoeffType cv_tr_over_R(unsigned int s) const;

       private:
        const ChemicalMixture<CoeffType> & _chem_mix;
        const ExternalThermo             & _ext_therm;
  };

  template <typename ExternalThermo, typename CoeffType>
  inline
  IdealGasInternalThermo<ExternalThermo, CoeffType>::IdealGasInternalThermo(const ExternalThermo & ext_therm, const ChemicalMixture<CoeffType> & chem_mix):
    _chem_mix(chem_mix),
    _ext_therm(ext_therm)
  {
     return;
  }

  template <typename ExternalThermo, typename CoeffType>
  inline
  IdealGasInternalThermo<ExternalThermo, CoeffType>::~IdealGasInternalThermo()
  {
     return;
  }

  template <typename ExternalThermo, typename CoeffType>
  inline
  const CoeffType IdealGasInternalThermo<ExternalThermo, CoeffType>::cv_tr_over_R(unsigned int s)
  {
     return _chem_mix.chemical_species[s]->n_tr_dofs();
  }

  template <typename ExternalThermo, typename CoeffType>
  inline
  const CoeffType IdealGasInternalThermo<ExternalThermo, CoeffType>::cv_tr(unsigned int s)
  {
    return _chem_mix.R(species) * (_chem_mix.chemical_species()[species])->n_tr_dofs();
  }

  template <typename ExternalThermo, typename CoeffType>
  inline
  const CoeffType IdealGasInternalThermo<ExternalThermo, CoeffType>::cv_rot(unsigned int s)
  {
    using std::max;

    return max(this->cv_tr(species) - this->cv_trans(species), CoeffType(0) ); 
  }

  template <typename ExternalThermo, typename CoeffType>
  inline
  const CoeffType IdealGasInternalThermo<ExternalThermo, CoeffType>::cv_rot_over_R(unsigned int s)
  {
    return std::max(this->cv_tr_over_R(s) - this->cv_trans_over_R(s), CoeffType(0) ); 
  }

  template <typename ExternalThermo, typename CoeffType>
  inline
  const CoeffType IdealGasInternalThermo<ExternalThermo, CoeffType>::cv_trans( const unsigned int species ) const
  {
    return CoeffType(1.5)*_chem_mixture.R(species);
  }

  template <typename ExternalThermo, typename CoeffType>
  inline
  const CoeffType IdealGasInternalThermo<ExternalThermo, CoeffType>::cv_trans_over_R( const unsigned int species ) const
  {
    return CoeffType(1.5);
  }


}


#endif
