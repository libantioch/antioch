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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef ANTIOCH_STAT_MECH_THERMO_H
#define ANTIOCH_STAT_MECH_THERMO_H

// Antioch
#include "antioch/chemical_mixture.h"
#include "antioch/input_utils.h"
#include "antioch/metaprogramming.h"
#include "antioch/antioch_exceptions.h"

// C++
#include <iomanip>
#include <vector>
#include <cmath>
#include <limits>

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
     * @returns species translational specific heat at constant volume.
     * Since the translational modes are assumed to be fully polulated
     * this is simply 
     * \f[
     *   C^{trans}_{v,s} \equiv \frac{\partial e^{trans}_s}{\partial T} = \frac{3}{2} R_s
     * \f]
     */
    CoeffType cv_trans( const unsigned int species ) const;

    /**
     * @returns species translational specific over R heat at constant volume.
     * Since the translational modes are assumed to be fully polulated
     * this is simply 
     * \f[
     *   \frac{C^{trans}_{v,s}}{\mathrm{R}} = \frac{3}{2}
     * \f]
     */
    CoeffType cv_trans_over_R( const unsigned int species ) const;


    /**
     * @returns species rotational specific heat at constant volume.
     * By convention, we lump the translational/rotational components
     * \f[
     *   C^{tr}_{v,s} \equiv C^{trans}_{v,s} + C^{rot}_{v,s}
     * \f]
     * so then
     * \f[
     *   C^{rot}_{v,s} \equiv C^{tr}_{v,s} - C^{trans}_{v,s}
     * \f]
     */
    CoeffType cv_rot( const unsigned int species ) const;

    /**
     * @returns species rotational specific heat at constant volume.
     * By convention, we lump the translational/rotational components
     * \f[
     *   C^{tr}_{v,s} \equiv C^{trans}_{v,s} + C^{rot}_{v,s}
     * \f]
     * so then
     * \f[
     *   \frac{C^{rot}_{v,s}}{\mathrm{R}} \equiv \frac{C^{tr}_{v,s}}{\mathrm{R}} - \frac{C^{trans}_{v,s}}{\mathrm{R}}
     * \f]
     */
    CoeffType cv_rot_over_R( const unsigned int species ) const;

    /**
     * @returns species translational/rotational specific heat at
     * constant volume.
     */
    CoeffType cv_tr (const unsigned int species) const;

    /**
     * @returns species translational/rotational specific heat over R at
     * constant volume.
     */
    CoeffType cv_tr_over_R (const unsigned int species) const;

      
    /**
     * @returns mixture translational/rotational specific heat at
     * constant volume.
     */
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type 
    >::type
    cv_tr (const VectorStateType& mass_fractions) const;

    /**
     * @returns species vibrational specific heat
     * constant volume.
     */
    template<typename StateType>
    StateType cv_vib (const unsigned int species, const StateType& Tv) const;

    /**
     * @returns species vibrational specific heat over R
     * constant volume.
     */
    template<typename StateType>
    StateType cv_vib_over_R (const unsigned int species, const StateType& Tv) const;
      
    /**
     * @returns mixture vibrational specific heat
     * constant volume.
     */
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type 
    >::type
    cv_vib (const typename Antioch::value_type<VectorStateType>::type& Tv, 
            const VectorStateType& mass_fractions) const;

    /**
     * @returns species electronic specific heat at constant volume.
     */
    template<typename StateType>
    StateType cv_el (const unsigned int species, const StateType& Te) const;
      
    /**
     * @returns mixture electronic specific heat at constant volume.
     */
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type 
    >::type
    cv_el (const typename Antioch::value_type<VectorStateType>::type& Te, 
           const VectorStateType& mass_fractions) const;

    /**
     * @returns species vibrational/electronic specific heat
     * constant volume.
     */
    template<typename StateType>
    StateType cv_ve (const unsigned int species, const StateType& Tv) const;
      
    /**
     * @returns mixture vibrational/electronic specific heat
     * constant volume.
     */
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type 
    >::type
    cv_ve (const typename Antioch::value_type<VectorStateType>::type& Tv, 
           const VectorStateType& mass_fractions) const;
      
    /**
     * @returns species total specific heat at constant volume, J/kg. Note that the input
     * translational/rotational temperature is currently not used, as these
     * modes are assumed to be fully excited.  However, the API is here
     * in case this assumption is removed later.
     */
    template<typename StateType>
    StateType cv (const unsigned int species, const StateType& T, const StateType& Tv) const;
      
    /**
     * @returns mixture total specific heat at constant volume, J/kg. Note that the input
     * translational/rotational temperature is currently not used, as these
     * modes are assumed to be fully excited.  However, the API is here
     * in case this assumption is removed later.
     */
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type 
    >::type
    cv (const typename Antioch::value_type<VectorStateType>::type& T, 
        const typename Antioch::value_type<VectorStateType>::type& Tv, 
        const VectorStateType& mass_fractions) const;

    /**
     * @returns mixture total specific heat at constant pressure, J/kg.
     */
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type 
    >::type
    cp (const typename Antioch::value_type<VectorStateType>::type& T, 
        const typename Antioch::value_type<VectorStateType>::type& Tv, 
        const VectorStateType& mass_fractions) const;

    /**
     * @returns species total enthalpy, J/kg.
     */
    template<typename StateType>
    StateType h_tot (const unsigned int species, const StateType& T, const StateType& Tv) const;

    /**
     * @returns species total enthalpy, J/kg.  Assumes thermal equilibrium
     * of translational/rotational and vibrational/electronic temperature.
     */
    template<typename StateType>
    StateType h_tot (const unsigned int species, const StateType& T) const;

    /**
     * @returns mixture specific enthalpy, J/kg.
     */
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type 
    >::type
    h_tot (const typename Antioch::value_type<VectorStateType>::type& T, 
           const typename Antioch::value_type<VectorStateType>::type& Tv, 
           const VectorStateType& mass_fractions) const;

    /**
     * @returns mixture specific enthalpy, J/kg. Assumes thermal equilibrium
     * of translational/rotational and vivbrational/electronic temperature.
     */
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type 
    >::type
    h_tot (const typename Antioch::value_type<VectorStateType>::type& T, 
           const VectorStateType& mass_fractions) const;
      
    /**
     * @returns species total internal energy, J/kg.
     */
    template<typename StateType>
    StateType e_tot (const unsigned int species, const StateType& T, const StateType& Tv) const;

    /**
     * @returns mixture total internal energy, J/kg.
     */
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type 
    >::type
    e_tot (const typename Antioch::value_type<VectorStateType>::type& T, 
           const typename Antioch::value_type<VectorStateType>::type& Tv, 
           const VectorStateType& mass_fractions) const;
      
    /**
     * @returns species total internal energy, J/kg.  Assumes thermal equilibrium
     * of translational/rotational and vibrational/electronic temperature.
     */
    template<typename StateType>
    StateType e_tot (const unsigned int species, const StateType& T) const;

    /**
     * @returns mixture total internal energy, J/kg.  Assumes thermal equilibrium
     * of translational/rotational and vivbrational/electronic temperature.
     */
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type 
    >::type
    e_tot (const typename Antioch::value_type<VectorStateType>::type& T, 
           const VectorStateType& mass_fractions) const;

    /**
     * @returns species translational/rotational energy, J/kg.
     */
    template<typename StateType>
    StateType e_tr (const unsigned int species, const StateType& T) const;
      
    /**
     * @returns mixture translational/rotational energy, J/kg.
     */
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type 
    >::type
    e_tr (const typename Antioch::value_type<VectorStateType>::type& T, 
          const VectorStateType& mass_fractions) const;

    /**
     * @returns species vibrational energy, J/kg.
     */
    template<typename StateType>
    StateType e_vib (const unsigned int species, const StateType& Tv) const;

    /**
     * @returns mixture vibrational energy, J/kg.
     */
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type 
    >::type
    e_vib (const typename Antioch::value_type<VectorStateType>::type& Tv, 
           const VectorStateType& mass_fractions) const;

    /**
     * @returns species electronic energy, J/kg.
     */
    template<typename StateType>
    StateType e_el (const unsigned int species, const StateType& Te) const;

    /**
     * @returns mixture electronic energy, J/kg.
     */
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type 
    >::type
    e_el (const typename Antioch::value_type<VectorStateType>::type& Te, 
          const VectorStateType& mass_fractions) const;
      
    /**
     * @returns species vibrational/electronic energy, J/kg.
     */
    template<typename StateType>
    StateType e_ve (const unsigned int species, const StateType& Tv) const;
      
    /**
     * @returns mixture vibrational/electronic energy, J/kg.
     */
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type 
    >::type
    e_ve (const typename Antioch::value_type<VectorStateType>::type& Te, 
          const VectorStateType& mass_fractions) const;

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
    template<typename VectorStateType>
    void e_and_cv_ve (const typename Antioch::value_type<VectorStateType>::type& Tv, 
                      const VectorStateType& mass_fractions,
                      typename Antioch::value_type<VectorStateType>::type &e_ve, 
                      typename Antioch::value_type<VectorStateType>::type &cv_ve) const;
      
    /**
     * @returns species formation energy, J/kg.
     */
    CoeffType e_0 (const unsigned int species) const;
      
    /**
     * @returns mixture formation energy, J/kg.
     */
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type 
    >::type
    e_0 (const VectorStateType& mass_fractions) const;
      
    /**
     * Computes the mixture vibrational temperature (K) from input
     * vibrational-electronic energy per unit mass (J/kg).
     */
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type 
    >::type
    Tv_from_e_ve (const typename Antioch::value_type<VectorStateType>::type& e_ve, 
                  const VectorStateType& mass_fractions,
                  typename Antioch::value_type<VectorStateType>::type Tv = -1) const;
      
    /**
     * Computes the mixture temperature (K) from input
     * total energy per unit mass (J/kg).
     */
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type 
    >::type
    T_from_e_tot (const typename Antioch::value_type<VectorStateType>::type& e_tot, 
                  const VectorStateType& mass_fractions,
                  typename Antioch::value_type<VectorStateType>::type T = -1) const;
      
    /**
     * Computes the mixture temperature (K) from input
     * translational/rotational energy per unit mass (J/kg).
     */
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type 
    >::type
    T_from_e_tr (const typename Antioch::value_type<VectorStateType>::type& e_tr, 
                 const VectorStateType& mass_fractions,
                 typename Antioch::value_type<VectorStateType>::type T = -1) const;
      
    /**
     * Computes the mixture temperature (K) from input
     * enthalpy per unit mass (J/kg).
     */
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type 
    >::type
    T_from_h_tot (const typename Antioch::value_type<VectorStateType>::type& h_tot, 
                  const VectorStateType& mass_fractions,
                  typename Antioch::value_type<VectorStateType>::type T = -1) const;
    
    /**
     * Computes the mixture temperature (K) from input
     * vibrational/electronic temperature Tv (K) and
     * enthalpy per unit mass (J/kg).
     */
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type 
    >::type
    T_from_h_tot_Tv (const typename Antioch::value_type<VectorStateType>::type& h_tot, 
                     const typename Antioch::value_type<VectorStateType>::type& Tv,
                     const VectorStateType& mass_fractions,
                     typename Antioch::value_type<VectorStateType>::type T = -1) const;
    
    /**
     * Computes the mixture specific entropy (J/kg-K) from input
     * temperature (K) and pressure (Pa).
     */
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type 
    >::type
    s (const typename Antioch::value_type<VectorStateType>::type& T, 
       const typename Antioch::value_type<VectorStateType>::type& p,
       const VectorStateType& mass_fractions) const;

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
    // NOP
  }

  template<typename CoeffType>
  inline
  StatMechThermodynamics<CoeffType>::~StatMechThermodynamics ()
  {
    // NOP
  }

  template<typename CoeffType>
  inline
  CoeffType StatMechThermodynamics<CoeffType>::cv_trans( const unsigned int species ) const
  {
    return CoeffType(1.5)*_chem_mixture.R(species);
  }

  template<typename CoeffType>
  inline
  CoeffType StatMechThermodynamics<CoeffType>::cv_trans_over_R( const unsigned int species ) const
  {
    return CoeffType(1.5);
  }

  template<typename CoeffType>
  inline
  CoeffType StatMechThermodynamics<CoeffType>::cv_rot( const unsigned int species ) const
  {
    using std::max;

    return max(this->cv_tr(species) - this->cv_trans(species), CoeffType(0) ); 
  }

  template<typename CoeffType>
  inline
  CoeffType StatMechThermodynamics<CoeffType>::cv_rot_over_R( const unsigned int species ) const
  {
    using std::max;

    return max(this->cv_tr_over_R(species) - this->cv_trans_over_R(species), CoeffType(0) ); 
  }

  template<typename CoeffType>
  inline
  CoeffType StatMechThermodynamics<CoeffType>::cv_tr (const unsigned int species) const
  {
    return _chem_mixture.R(species)*(_chem_mixture.chemical_species()[species])->n_tr_dofs();
  }

  template<typename CoeffType>
  inline
  CoeffType StatMechThermodynamics<CoeffType>::cv_tr_over_R (const unsigned int species) const
  {
    return (_chem_mixture.chemical_species()[species])->n_tr_dofs();
  }

  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type 
  >::type
  StatMechThermodynamics<CoeffType>::cv_tr (const VectorStateType& mass_fractions) const
  {
    typename Antioch::value_type<VectorStateType>::type 
      cv_tr = mass_fractions[0]*this->cv_tr(0);

    for( unsigned int s = 1; s < _chem_mixture.n_species(); s++ )
      {
        cv_tr += mass_fractions[s]*this->cv_tr(s);
      }
    
    return cv_tr;
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::cv_vib (const unsigned int species, 
                                                       const StateType& Tv) const
  {
      return this->cv_vib_over_R(species,Tv) * (_chem_mixture.chemical_species()[species])->gas_constant();
  }
      
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::cv_vib_over_R (const unsigned int species, 
                                                       const StateType& Tv) const
  {
    using std::exp;

    // convenience
    const ChemicalSpecies<CoeffType>& chem_species = *(_chem_mixture.chemical_species()[species]);
    const std::vector<CoeffType>& theta_v  = chem_species.theta_v();
    const std::vector<unsigned int>& ndg_v = chem_species.ndg_v();
    
    antioch_assert_equal_to(ndg_v.size(), theta_v.size());
    
    // Use an input datum to make sure we get the size right
    StateType cv_vib = Antioch::zero_clone(Tv);
	
    if (theta_v.empty())
      return cv_vib;

    for (unsigned int level=0; level<ndg_v.size(); level++)
      {
        typedef typename Antioch::raw_value_type<StateType>::type raw_type;
        const StateType expval = exp(theta_v[level]/Tv);
        const StateType expvalminusone = expval - raw_type(1);
      
        cv_vib += (static_cast<CoeffType>(ndg_v[level])*
                   theta_v[level]*theta_v[level]*expval/(expvalminusone*expvalminusone)/(Tv*Tv));
      }
    
    return cv_vib;
  }
      
  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type 
  >::type
  StatMechThermodynamics<CoeffType>::cv_vib (const typename Antioch::value_type<VectorStateType>::type& Tv, 
                                             const VectorStateType& mass_fractions) const
  {
    typename Antioch::value_type<VectorStateType>::type
      cv_vib = mass_fractions[0]*this->cv_vib(0, Tv);

    for( unsigned int s = 1; s < _chem_mixture.n_species(); s++ )
      {
        cv_vib += mass_fractions[s]*this->cv_vib(s, Tv);
      }
    
    return cv_vib;
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::cv_el (const unsigned int species, 
                                                      const StateType& Te) const
  {
    using std::exp;

    // convenience
    const ChemicalSpecies<CoeffType>& chem_species = *(_chem_mixture.chemical_species()[species]);
    const std::vector<CoeffType>& theta_e  = chem_species.theta_e();
    const std::vector<unsigned int>& ndg_e = chem_species.ndg_e();
    
    antioch_assert_equal_to(ndg_e.size(), theta_e.size());
    
    StateType cv_el = Antioch::zero_clone(Te);
    
    // Really < 2?  Yes, b/c theta_e[0] = 0.0 always.  See
    // read_species_electronic_data_ascii_default in
    // species_ascii_parsing.h
    if (theta_e.size() < 2) return cv_el;
    
    typedef typename Antioch::raw_value_type<StateType>::type raw_type;
    const raw_type one = static_cast<raw_type>(1);

    const StateType Teinv = one/Te;
    const StateType Te2inv = Teinv*Teinv;
    
    StateType
      num = Antioch::zero_clone(Te), dnum = Antioch::zero_clone(Te),
      den = Antioch::zero_clone(Te), dden = Antioch::zero_clone(Te);
    
    for (unsigned int level=0; level<theta_e.size(); level++)
      {
        const StateType 
          expval = exp (-theta_e[level] * Teinv),
          den_l  = static_cast<raw_type>(ndg_e[level])*expval,
          num_l  = den_l*theta_e[level],
          dden_l = num_l*Te2inv,
          dnum_l = dden_l*theta_e[level];
        
        num  += num_l;
        den  += den_l;
	
        dden += dden_l;
        dnum += dnum_l;
      }
    
    const StateType invden = one/den;
    
    cv_el = chem_species.gas_constant()*(dnum - num*dden*invden) * invden;

    return cv_el;
  }
  
  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type 
  >::type
  StatMechThermodynamics<CoeffType>::cv_el (const typename Antioch::value_type<VectorStateType>::type& Te, 
                                            const VectorStateType& mass_fractions) const
  {
    typename Antioch::value_type<VectorStateType>::type
      cv_el = mass_fractions[0]*this->cv_el(0, Te);

    for( unsigned int s = 1; s < _chem_mixture.n_species(); s++ )
      {
        cv_el += mass_fractions[s]*this->cv_el(s, Te);
      }
    
    return cv_el;
  }
      
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::cv_ve (const unsigned int species, 
                                                      const StateType& Tv) const
  {
    return (this->cv_vib(species, Tv) + this->cv_el(species, Tv));
  }
      
  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type 
  >::type
  StatMechThermodynamics<CoeffType>::cv_ve (const typename Antioch::value_type<VectorStateType>::type& Tv, 
                                            const VectorStateType& mass_fractions) const

  {
    return (this->cv_vib(Tv, mass_fractions) + this->cv_el(Tv, mass_fractions));
  }
      
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::cv (const unsigned int species, 
                                                   const StateType& /* T */, 
                                                   const StateType& Tv) const
  {
    return (this->cv_tr(species) + this->cv_ve(species, Tv));
  }
      
  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type 
  >::type
  StatMechThermodynamics<CoeffType>::cv (const typename Antioch::value_type<VectorStateType>::type& /* T */, 
                                         const typename Antioch::value_type<VectorStateType>::type& Tv, 
                                         const VectorStateType& mass_fractions) const
  {
    return (this->cv_tr(mass_fractions) + this->cv_ve(Tv, mass_fractions));
  }

  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type 
  >::type
  StatMechThermodynamics<CoeffType>::cp (const typename Antioch::value_type<VectorStateType>::type& T, 
                                         const typename Antioch::value_type<VectorStateType>::type& Tv, 
                                         const VectorStateType& mass_fractions) const
  {
    return (this->cv(T,Tv,mass_fractions) + this->_chem_mixture.R(mass_fractions));
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::h_tot (const unsigned int species, 
                                                      const StateType& T, 
                                                      const StateType& Tv) const
  {
    const ChemicalSpecies<CoeffType>& chem_species = *(_chem_mixture.chemical_species()[species]);
    return (this->e_tot(species, T, Tv) + chem_species.gas_constant()*T);
  }


  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::h_tot (const unsigned int species, 
                                                      const StateType& T) const
  {
    return this->h_tot(species, T, T);
  }

  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type 
  >::type
  StatMechThermodynamics<CoeffType>::h_tot (const typename Antioch::value_type<VectorStateType>::type& T, 
                                            const typename Antioch::value_type<VectorStateType>::type& Tv, 
                                            const VectorStateType& mass_fractions) const
  {
    typename Antioch::value_type<VectorStateType>::type
      h_tot = mass_fractions[0]*this->h_tot(0, T, Tv);

    for( unsigned int s = 1; s < _chem_mixture.n_species(); s++ )
      {
        h_tot += mass_fractions[s]*this->h_tot(s, T, Tv);
      }
    
    return h_tot;
  }

  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type 
  >::type
  StatMechThermodynamics<CoeffType>::h_tot (const typename Antioch::value_type<VectorStateType>::type& T, 
                                            const VectorStateType& mass_fractions) const
  {
    return this->h_tot(T, T, mass_fractions);
  }
      
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::e_tot (const unsigned int species, 
                                                      const StateType& T, 
                                                      const StateType& Tv) const
  {
    return (this->e_tr(species, T) + this->e_ve(species, Tv) + this->e_0(species));
  }


  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type 
  >::type
  StatMechThermodynamics<CoeffType>::e_tot (const typename Antioch::value_type<VectorStateType>::type& T, 
                                            const typename Antioch::value_type<VectorStateType>::type& Tv, 
                                            const VectorStateType& mass_fractions) const
  {
    return (this->e_tr(T, mass_fractions) + this->e_ve(Tv, mass_fractions) + this->e_0(mass_fractions));
  }
      
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::e_tot (const unsigned int species, 
                                                      const StateType& T) const
  {
    return this->e_tot(species, T, T);
  }

  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type 
  >::type
  StatMechThermodynamics<CoeffType>::e_tot (const typename Antioch::value_type<VectorStateType>::type& T, 
                                            const VectorStateType& mass_fractions) const
  {
    return this->e_tot(T, T, mass_fractions);
  }
      
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::e_tr (const unsigned int species, 
                                                     const StateType& T) const
  {
    return this->cv_tr(species)*T;
  }
      
  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type 
  >::type
  StatMechThermodynamics<CoeffType>::e_tr (const typename Antioch::value_type<VectorStateType>::type& T, 
                                           const VectorStateType& mass_fractions) const
  {
    typename Antioch::value_type<VectorStateType>::type
      e_tr = mass_fractions[0]*this->e_tr(0, T);

    for( unsigned int s = 1; s < _chem_mixture.n_species(); s++ )
      {
        e_tr += mass_fractions[s]*this->e_tr(s, T);
      }
    
    return e_tr;
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::e_vib (const unsigned int species, 
                                                      const StateType& Tv) const
  {
    using std::exp;

    // convenience
    const ChemicalSpecies<CoeffType>& chem_species = *(_chem_mixture.chemical_species()[species]);
    const std::vector<CoeffType>& theta_v  = chem_species.theta_v();
    const std::vector<unsigned int>& ndg_v = chem_species.ndg_v();
    
    antioch_assert_equal_to(ndg_v.size(), theta_v.size());
    
    StateType e_vib = 0.0;

    if (theta_v.empty()) return e_vib;
    
    for (unsigned int level=0; level<ndg_v.size(); level++)
      e_vib += ndg_v[level]*chem_species.gas_constant()*theta_v[level]/(exp(theta_v[level]/Tv) - 1.);
    
    return e_vib;
  }

  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type 
  >::type
  StatMechThermodynamics<CoeffType>::e_vib (const typename Antioch::value_type<VectorStateType>::type& Tv, 
                                            const VectorStateType& mass_fractions) const
  {
    typename Antioch::value_type<VectorStateType>::type
      e_vib = mass_fractions[0]*this->e_vib(0, Tv);

    for( unsigned int s = 1; s < _chem_mixture.n_species(); s++ )
      {
        e_vib += mass_fractions[s]*this->e_vib(s, Tv);
      }
    
    return e_vib;
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::e_el (const unsigned int species,
                                                     const StateType& Te) const
  {
    using std::exp;

    // convenience
    const ChemicalSpecies<CoeffType>& chem_species = *(_chem_mixture.chemical_species()[species]);
    const std::vector<CoeffType>& theta_e  = chem_species.theta_e();
    const std::vector<unsigned int>& ndg_e = chem_species.ndg_e();
    
    antioch_assert_equal_to(ndg_e.size(), theta_e.size());

    StateType e_el = Antioch::zero_clone(Te);

    if (theta_e.size() < 2) return e_el;

    StateType num = Antioch::zero_clone(Te), den = Antioch::zero_clone(Te);
    
    for (unsigned int level=0; level<theta_e.size(); level++)
      {
        const StateType expval = exp (-theta_e[level] / Te);
        num += static_cast<StateType>(ndg_e[level])*theta_e[level]*expval;
        den += static_cast<StateType>(ndg_e[level])*expval;	  
      }
    
    return chem_species.gas_constant() * num / den;
  }

  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type 
  >::type
  StatMechThermodynamics<CoeffType>::e_el (const typename Antioch::value_type<VectorStateType>::type& Te, 
                                           const VectorStateType& mass_fractions) const
  {
    typename Antioch::value_type<VectorStateType>::type
      e_el = mass_fractions[0]*this->e_el(0, Te);

    for( unsigned int s = 1; s < _chem_mixture.n_species(); s++ )
      {
        e_el += mass_fractions[s]*this->e_el(s, Te);
      }
    
    return e_el;
  }
      
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType StatMechThermodynamics<CoeffType>::e_ve (const unsigned int species, 
                                                     const StateType& Tv) const
  {
    return (this->e_vib(species,Tv) + this->e_el(species,Tv));
  }
      
  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type 
  >::type
  StatMechThermodynamics<CoeffType>::e_ve (const typename Antioch::value_type<VectorStateType>::type& Tv, 
                                           const VectorStateType& mass_fractions) const
  {
    return (this->e_vib(Tv,mass_fractions) + this->e_el(Tv,mass_fractions));
  }

  template<typename CoeffType>
  inline
  CoeffType StatMechThermodynamics<CoeffType>::e_ve_min () const
  {
    antioch_not_implemented();
  }

  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  void StatMechThermodynamics<CoeffType>::e_and_cv_ve (const typename Antioch::value_type<VectorStateType>::type& Tv, 
                                                       const VectorStateType& mass_fractions,
                                                       typename Antioch::value_type<VectorStateType>::type &e_ve, 
                                                       typename Antioch::value_type<VectorStateType>::type &cv_ve) const
  {
    e_ve  = this->e_ve (Tv, mass_fractions);
    cv_ve = this->cv_ve(Tv, mass_fractions);
  }
      
  template<typename CoeffType>
  inline
  CoeffType StatMechThermodynamics<CoeffType>::e_0 (const unsigned int species) const
  {
    const ChemicalSpecies<CoeffType>& chem_species = *(_chem_mixture.chemical_species()[species]);
    return chem_species.formation_enthalpy();
  }
      
  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type 
  >::type
  StatMechThermodynamics<CoeffType>::e_0 (const VectorStateType& mass_fractions) const
  {
    typename Antioch::value_type<VectorStateType>::type
      e_0 = mass_fractions[0]*this->e_0(0);
    
    for( unsigned int s = 1; s < _chem_mixture.n_species(); s++ )
      {
        e_0 += mass_fractions[s]*this->e_0(s);
      }
    
    return e_0;
  }
      
  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type 
  >::type
  StatMechThermodynamics<CoeffType>::Tv_from_e_ve 
    (const typename Antioch::value_type<VectorStateType>::type& e_ve, 
     const VectorStateType& mass_fractions,
     typename Antioch::value_type<VectorStateType>::type Tv) const
  {
    antioch_not_implemented();
  }
      
  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type 
  >::type
  StatMechThermodynamics<CoeffType>::T_from_e_tot (const typename Antioch::value_type<VectorStateType>::type& e_tot, 
                                                   const VectorStateType& mass_fractions,
                                                   typename Antioch::value_type<VectorStateType>::type T) const
  {
    using std::abs;
    using std::max;
    using std::min;

    typedef typename Antioch::value_type<VectorStateType>::type StateType;

    // Cache the translational/rotational specific heat - this will be used repeatedly 
    // and involves (2*NS-1) flops to compute, and since this has no functional 
    // dependence on temperature it will not change throughout the Newton iteration.
    const typename Antioch::value_type<VectorStateType>::type
      Cv_tr = this->cv_tr(mass_fractions);

    // Similarly for mixture formation energy
    const typename Antioch::value_type<VectorStateType>::type
      E_0 = this->e_0(mass_fractions);

    // if the user does not provide an initial guess for the temperature
    // assume it is all in translation/rotation to compute a starting value.
    if (T < 0)
      {
	T = (e_tot - E_0) / Cv_tr;
	T = min(max(T,StateType(10.)),StateType(20000.));
	
        // FIXME: Use Antioch::Limits or similar? (i.e., don't
        // hardcode min and max T)

	// make sure the initial guess is valid
	//T = max(T, Limits::CompNSLimits::T_min());
        T = max(T, StateType(10.));
	T = min(T, StateType(2.e4));
      }
    
    // compute the translational/rotational temperature of the mixture using Newton-Rhapson iteration
    CoeffType delta_T = std::numeric_limits<StateType>::max();
    const unsigned int max_iterations = 100;

    const CoeffType dT_reltol= std::numeric_limits<CoeffType>::epsilon() * 100;
    
    // NOTE: FIN-S uses a hardcoded, absolute tolerance on delta_T of
    // 1e-8.  Using a relative tolerance here of 100*epsilon.
    for (unsigned int iter = 0;
         abs(delta_T/T) > dT_reltol &&
         T >= 0.;
         ++iter)
      {
	if (iter == max_iterations)
	  throw FailedNewtonTTvInversion ("ERROR: failed to converge T_from_e_tot!");
	
	// compute the residual, defined as the mismatch between the input e_tot and
	// the value corresponding to the current T iterate
        typename Antioch::value_type<VectorStateType>::type
          Re = 0., dRevdT = 0.;
        
        this->e_and_cv_ve(T, mass_fractions, Re, dRevdT);
        Re += this->e_tr(T, mass_fractions) + E_0 - e_tot;
        dRevdT += Cv_tr;

	delta_T = -Re / dRevdT;
	T += delta_T;
      }

    return T;
  }
      
  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type 
  >::type
  StatMechThermodynamics<CoeffType>::T_from_e_tr
    (const typename Antioch::value_type<VectorStateType>::type& e_tr, 
     const VectorStateType& mass_fractions,
     typename Antioch::value_type<VectorStateType>::type T ) const
  {
    antioch_not_implemented();
  }
      
  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type 
  >::type
  StatMechThermodynamics<CoeffType>::T_from_h_tot
    (const typename Antioch::value_type<VectorStateType>::type& h_tot, 
     const VectorStateType& mass_fractions,
     typename Antioch::value_type<VectorStateType>::type T) const
  {
    antioch_not_implemented();
  }
    
  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type 
  >::type
  StatMechThermodynamics<CoeffType>::T_from_h_tot_Tv
    (const typename Antioch::value_type<VectorStateType>::type& h_tot, 
     const typename Antioch::value_type<VectorStateType>::type& Tv,
     const VectorStateType& mass_fractions,
     typename Antioch::value_type<VectorStateType>::type T) const
  {
    antioch_not_implemented();
  }

    
  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type 
  >::type
  StatMechThermodynamics<CoeffType>::s
    (const typename Antioch::value_type<VectorStateType>::type& T, 
     const typename Antioch::value_type<VectorStateType>::type& p,
     const VectorStateType& mass_fractions) const
  {
    antioch_not_implemented();
  }

} // end namespace Antioch

#endif // ANTIOCH_STAT_MECH_THERMO_H
