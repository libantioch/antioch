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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef ANTIOCH_STAT_MECH_THERMO_H
#define ANTIOCH_STAT_MECH_THERMO_H

// C++
#include <iomanip>
#include <vector>
#include <cmath>

// Antioch
#include "antioch/chemical_mixture.h"
#include "antioch/input_utils.h"
#include "antioch/metaprogramming.h"

namespace Antioch
{

  template<typename CoeffType=double>
  class StatMechThermodynamics
  {
  public:
    
    StatMechThermodynamics( const ChemicalMixture<CoeffType>& chem_mixture );

    //! Destructor
    /*! virtual so this can be subclassed by the user. */
    virtual ~StatMechThermodynamics();

    /**
     * @returns species translational/rotational specific heat at
     * constant volume.
     */
    CoeffType cv_tr (const unsigned int species) const;

      
    /**
     * @returns mixture translational/rotational specific heat at
     * constant volume.
     */
    CoeffType cv_tr () const;
      
    /**
     * @returns species vibrational specific heat
     * constant volume.
     */
    template<typename StateType>
    StateType cv_vib (const unsigned int species, const StateType Tv) const;
      
    /**
     * @returns mixture vibrational specific heat
     * constant volume.
     */
    template<typename StateType>
    StateType cv_vib (const StateType Tv) const;
      
    /**
     * @returns species vibrational/electronic specific heat
     * constant volume.
     */
    template<typename StateType>
    StateType cv_ve (const unsigned int species, const StateType Tv) const;
      
    /**
     * @returns mixture vibrational/electronic specific heat
     * constant volume.
     */
    template<typename StateType>
    StateType cv_ve (const StateType Tv) const;
      
    /**
     * @returns species total specific heat at constant volume, J/kg. Note that the input
     * translational/rotational temperature is currently not used, as these
     * modes are assumed to be fully excited.  However, the API is here
     * in case this assumption is removed later.
     */
    template<typename StateType>
    StateType cv (const unsigned int species, const StateType T, const StateType Tv) const;
      
    /**
     * @returns mixture total specific heat at constant volume, J/kg. Note that the input
     * translational/rotational temperature is currently not used, as these
     * modes are assumed to be fully excited.  However, the API is here
     * in case this assumption is removed later.
     */
    template<typename StateType>
    StateType cv (const StateType T, const StateType Tv) const;

    /**
     * @returns mixture total specific heat at constant pressure, J/kg.
     */
    template<typename StateType>
    StateType cp (const StateType Ttr, const StateType Tv) const;

    /**
     * @returns species total enthalpy, J/kg.
     */
    template<typename StateType>
    StateType h_tot (const unsigned int species, const StateType T, const StateType Tv) const;

    /**
     * @returns species total enthalpy, J/kg.  Assumes thermal equilibrium
     * of translational/rotational and vibrational/electronic temperature.
     */
    template<typename StateType>
    StateType h_tot (const unsigned int species, const StateType T) const;

    /**
     * @returns mixture specific enthalpy, J/kg.
     */
    template<typename StateType>
    StateType h_tot (const StateType T, const StateType Tv) const;

    /**
     * @returns mixture specific enthalpy, J/kg. Assumes thermal equilibrium
     * of translational/rotational and vivbrational/electronic temperature.
     */
    template<typename StateType>
    StateType h_tot (const StateType T) const;
      
    /**
     * @returns species total internal energy, J/kg.
     */
    template<typename StateType>
    StateType e_tot (const unsigned int species, const StateType T, const StateType Tv) const;

    /**
     * @returns mixture total internal energy, J/kg.
     */
    template<typename StateType>
    StateType e_tot (const StateType T, const StateType Tv) const;
      
    /**
     * @returns species total internal energy, J/kg.  Assumes thermal equilibrium
     * of translational/rotational and vibrational/electronic temperature.
     */
    template<typename StateType>
    StateType e_tot (const unsigned int species, const StateType T) const;

    /**
     * @returns mixture total internal energy, J/kg.  Assumes thermal equilibrium
     * of translational/rotational and vivbrational/electronic temperature.
     */
    template<typename StateType>
    StateType e_tot (const StateType T) const;
      
    /**
     * @returns species translational/rotational energy, J/kg.
     */
    template<typename StateType>
    StateType e_tr (const unsigned int species, const StateType T) const;
      
    /**
     * @returns mixture translational/rotational energy, J/kg.
     */
    template<typename StateType>
    StateType e_tr (const StateType T) const;

    /**
     * @returns species vibrational/electronic energy, J/kg.
     */
    template<typename StateType>
    StateType e_ve (const unsigned int species, const StateType Tv) const;
      
    /**
     * @returns mixture vibrational energy, J/kg.
     */
    template<typename StateType>
    StateType e_vib (const StateType Tv) const;
      
    /**
     * @returns species vibrational energy, J/kg.
     */
    template<typename StateType>
    StateType e_vib (const unsigned int species, const StateType Tv) const;
      
    /**
     * @returns mixture vibrational/electronic energy, J/kg.
     */
    template<typename StateType>
    StateType e_ve (const StateType Tv) const;

    /**
     * @returns the minimum valid mixture vibrational/electronic energy, J/kg.
     * We allow for a user-specified minimum vibrational temperature, Tv.
     * This in turn sets an a priori minimum valid species 
     * vibrational/electronic temperature.  However, the minimum valid
     * value for a mixture depends on the species mass fractions and must
     * be computed.  That is what this method does.
     */
    CoeffType e_ve_min () const;

    /**
     * @returns mixture vibrational/electronic energy (same as
     * e_ve()) and specific heat (same as cv_ve()), calculated
     * together for efficiency.
     */
    template<typename StateType>
    void e_and_cv_ve (const StateType Tv, StateType &e_ve, StateType &cv_ve) const;
      
    /**
     * @returns species formation energy, J/kg.
     */
    CoeffType e_0 (const unsigned int species) const;
      
    /**
     * @returns mixture formation energy, J/kg.
     */
    CoeffType e_0 () const;
      
    /**
     * Computes the mixture vibrational temperature (K) from input
     * vibrational-electronic energy per unit mass (J/kg).
     */
    template<typename StateType>
    StateType Tv_from_e_ve (const StateType e_ve, StateType Tv = -1.) const;
      
    /**
     * Computes the mixture temperature (K) from input
     * total energy per unit mass (J/kg).
     */
    template<typename StateType>
    StateType T_from_e_tot (const StateType e_tot, StateType T = -1.) const;
      
    /**
     * Computes the mixture temperature (K) from input
     * translational/rotational energy per unit mass (J/kg).
     */
    template<typename StateType>
    StateType T_from_e_tr (const StateType e_tr, StateType T = -1.) const;
      
    /**
     * Computes the mixture temperature (K) from input
     * enthalpy per unit mass (J/kg).
     */
    template<typename StateType>
    StateType T_from_h_tot (const StateType h_tot, StateType T = -1.) const;
    
    /**
     * Computes the mixture temperature (K) from input
     * vibrational/electronic temperature Tv (K) and
     * enthalpy per unit mass (J/kg).
     */
    template<typename StateType>
    StateType T_from_h_tot_Tv (const StateType h_tot, const StateType Tv, StateType T = -1.) const;
    
    /**
     * Computes the mixture specific entropy (J/kg-K) from input
     * temperature (K) and pressure (Pa).
     */
    template<typename StateType>
    StateType s (const StateType T, const StateType p ) const;

  protected:

    const ChemicalMixture<CoeffType>& _chem_mixture;

  private:
    
    //! Default constructor
    /*! Private to force to user to supply a ChemicalMixture object.*/
    StatMechThermodynamics();
  };

  
  /* ------------------------- Inline Functions -------------------------*/
  
  template<typename CoeffType>
  inline
  StatMechThermodynamics<CoeffType>::StatMechThermodynamics( const ChemicalMixture<CoeffType>& chem_mixture )
    : _chem_mixture(chem_mixture)
  {
    /*! \todo Read all necessary data */
  }

  template<typename CoeffType>
  inline
  StatMechThermodynamics<CoeffType>::~StatMechThermodynamics ()
  {
    // NOP
  }

  template<typename CoeffType>
  inline
  CoeffType StatMechThermodynamics<CoeffType>::cv_tr (const unsigned int species) const
  {
    return _chem_mixture.R(species)*(_chem_mixture.chemical_species()[species])->n_tr_dofs();
  }

  template<typename CoeffType>
  inline
  CoeffType StatMechThermodynamics<CoeffType>::cv_tr () const
  {
    antioch_not_implemented();
  }
      
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::cv_vib (const unsigned int species, 
                                                       const StateType Tv) const
  {
    antioch_not_implemented();
  }
      
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::cv_vib (const StateType Tv) const
  {
    antioch_not_implemented();
  }
      
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::cv_ve (const unsigned int species, 
                                                      const StateType Tv) const
  {
    antioch_not_implemented();
  }
      
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::cv_ve (const StateType Tv) const
  {
    antioch_not_implemented();
  }
      
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::cv (const unsigned int species, 
                                                   const StateType T, 
                                                   const StateType Tv) const
  {
    antioch_not_implemented();
  }
      
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::cv (const StateType T, 
                                                   const StateType Tv) const
  {
    antioch_not_implemented();
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::cp (const StateType Ttr, 
                                                   const StateType Tv) const
  {
    antioch_not_implemented();
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::h_tot (const unsigned int species, 
                                                      const StateType T, 
                                                      const StateType Tv) const
  {
    antioch_not_implemented();
  }


  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::h_tot (const unsigned int species, 
                                                      const StateType T) const
  {
    antioch_not_implemented();
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::h_tot (const StateType T, 
                                                      const StateType Tv) const
  {
    antioch_not_implemented();
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::h_tot (const StateType T) const
  {
    antioch_not_implemented();
  }
      
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::e_tot (const unsigned int species, 
                                                      const StateType T, 
                                                      const StateType Tv) const
  {
    antioch_not_implemented();
  }


  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::e_tot (const StateType T, 
                                                      const StateType Tv) const
  {
    antioch_not_implemented();
  }
      
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::e_tot (const unsigned int species, 
                                                      const StateType T) const
  {
    antioch_not_implemented();
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::e_tot (const StateType T) const
  {
    antioch_not_implemented();
  }
      
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::e_tr (const unsigned int species, 
                                                     const StateType T) const
  {
    antioch_not_implemented();
  }
      
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::e_tr (const StateType T) const
  {
    antioch_not_implemented();
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::e_ve (const unsigned int species, 
                                                     const StateType Tv) const
  {
    antioch_not_implemented();
  }
      
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::e_vib (const StateType Tv) const
  {
    antioch_not_implemented();
  }
      
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::e_vib (const unsigned int species, 
                                                      const StateType Tv) const
  {
    antioch_not_implemented();
  }
      
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::e_ve (const StateType Tv) const
  {
    antioch_not_implemented();
  }

  template<typename CoeffType>
  inline
  CoeffType StatMechThermodynamics<CoeffType>::e_ve_min () const
  {
    antioch_not_implemented();
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  void StatMechThermodynamics<CoeffType>::e_and_cv_ve (const StateType Tv, 
                                                       StateType &e_ve, 
                                                       StateType &cv_ve) const
  {
    antioch_not_implemented();
  }
      
  template<typename CoeffType>
  inline
  CoeffType StatMechThermodynamics<CoeffType>::e_0 (const unsigned int species) const
  {
    antioch_not_implemented();
  }
      
  template<typename CoeffType>
  inline
  CoeffType StatMechThermodynamics<CoeffType>::e_0 () const
  {
    antioch_not_implemented();
  }
      
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::Tv_from_e_ve (const StateType e_ve, 
                                                             StateType Tv ) const
  {
    antioch_not_implemented();
  }
      
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::T_from_e_tot (const StateType e_tot, 
                                                             StateType T ) const
  {
    antioch_not_implemented();
  }
      
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::T_from_e_tr (const StateType e_tr, 
                                                            StateType T ) const
  {
    antioch_not_implemented();
  }
      
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::T_from_h_tot (const StateType h_tot, 
                                                             StateType T ) const
  {
    antioch_not_implemented();
  }
    
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::T_from_h_tot_Tv (const StateType h_tot, 
                                                                const StateType Tv, 
                                                                StateType T ) const
  {
    antioch_not_implemented();
  }

    
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::s (const StateType T, 
                                                  const StateType p ) const
  {
    antioch_not_implemented();
  }

} // end namespace Antioch

#endif // ANTIOCH_STAT_MECH_THERMO_H
