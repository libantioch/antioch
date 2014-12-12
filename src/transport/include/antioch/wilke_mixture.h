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

#ifndef ANTIOCH_WILKE_MIXTURE_H
#define ANTIOCH_WILKE_MIXTURE_H

// Antioch
#include "antioch/physics_metaprogramming_decl.h"
#include "antioch/temp_cache.h"
#include "antioch/cmath_shims.h"

// C++
#include <vector>

namespace Antioch
{
  /*
      This somewhat ridicule thermodynamics
      construct is because StatMechThermodynamics
      can't give cp, and NASA9 can't give
      the vibrational, translational, etc.,
      parts of this same quantity.
      Notice that the internal R used is the
      massic one, thus Cp is in J/kg/K

      \todo make Wilke rules be computed here?
  */
  template<class Diffusion, class Viscosity, class ThermalConduction, // physics
           class Mixture, class ThermoEvaluator,                      // mixture + thermo
           class CoeffType = double>                                  // type
  class WilkeMixture
  {
  public:

    WilkeMixture( const Diffusion & diffusion, const Viscosity & viscosity, 
                  const ThermalConduction & thermal_conduction,
                  const Mixture & transport_mixture,
                  const ThermoEvaluator & thermo_eval );
    ~WilkeMixture();

    //! const ref to transport_mixture
    const Mixture & transport_mixture() const; // contains the detailed thermo for species

    //! const ref to thermo evaluator, internal computations only (vib, rot, trans)
    const ThermoEvaluator & thermo_evaluator() const; // contains the non detailed thermo for mixture

    //! const ref to diffusion
    const Diffusion & diffusion() const;

    //! const ref to viscosity
    const Viscosity & viscosity() const;

    //! const ref to thermal conduction
    const ThermalConduction & thermal_conduction() const;


    //! Full array of mixture-level diffusion
    template <typename StateType, typename VectorStateType>
    void D(const StateType & T, const StateType & rho, const StateType & cTot, 
           const VectorStateType & mass_fractions, const VectorStateType & mu, VectorStateType & ds) const;

    //! Viscosity of one species
    template <typename StateType>
    const ANTIOCH_AUTO(StateType)
    mu(unsigned int s, const StateType & T) const;

    //! Viscosities of all species
    template <typename StateType, typename VectorStateType>
    void mu(const StateType & T, VectorStateType & mu) const;

    //! Thermal conduction of one species
    template <typename StateType>
    const ANTIOCH_AUTO(StateType)
    k(unsigned int s, const StateType & mu, const StateType & T, const StateType &n_molar_mixture = 0) const;

    //! Thermal conduction of all species
    template <typename StateType, typename VectorStateType>
    void k(const VectorStateType & mu, const StateType & T, const VectorStateType & mass_fractions, const StateType &rho, VectorStateType & k) const;

    //! Thermal conduction and diffusion
    template <typename StateType, typename VectorStateType>
    void D_and_k(const VectorStateType & mu, const StateType & T, const StateType & rho, 
                 const VectorStateType & mass_fractions, VectorStateType & k, VectorStateType & D) const;


    CoeffType Mr_Ms_to_the_one_fourth( const unsigned int r,
                                       const unsigned int s ) const;
    
    CoeffType denominator( const unsigned int r,
                           const unsigned int s ) const;

  protected:

    const Mixture           & _transport_mixture;

    const Diffusion         & _diffusion_set;

    const Viscosity         & _viscosity_set;

    const ThermalConduction & _thermal_conduction_set;

    const ThermoEvaluator   & _thermo_evaluator;

    //! Cache for numerator term
    /*! \todo We should use a more efficient data structure */
    std::vector<std::vector<CoeffType> > _Mr_Ms_to_the_one_fourth;
    
    //! Cache for denominator term
    /*! \todo We should use a more efficient data structure */
    std::vector<std::vector<CoeffType> > _denom;

  };

  template<typename Diffusion, typename Viscosity, typename ThermalConduction, typename Mixture, typename ThermoEvaluator, typename CoeffType>
  inline
  WilkeMixture<Diffusion,Viscosity,ThermalConduction,Mixture,ThermoEvaluator,CoeffType>::WilkeMixture(
                                        const Diffusion & diffusion, const Viscosity & viscosity, 
                                        const ThermalConduction & thermal_conduction,
                                        const Mixture & transport_mixture,
                                        const ThermoEvaluator & thermo_eval )
    : _transport_mixture(transport_mixture),
      _diffusion_set(diffusion),
      _viscosity_set(viscosity),
      _thermal_conduction_set(thermal_conduction),
      _thermo_evaluator(thermo_eval),
      _Mr_Ms_to_the_one_fourth(_transport_mixture.n_species(), std::vector<CoeffType>(_transport_mixture.n_species())),
      _denom(_transport_mixture.n_species(), std::vector<CoeffType>(_transport_mixture.n_species()))
  {

    for( unsigned int r = 0; r < transport_mixture.n_species(); r++ )
      {
        for( unsigned int s = 0; s < transport_mixture.n_species(); s++ )
          {
            const CoeffType Mr = _transport_mixture.chemical_mixture().M(r);
            const CoeffType Ms = _transport_mixture.chemical_mixture().M(s);

            _Mr_Ms_to_the_one_fourth[r][s] = ant_pow( Mr/Ms, Antioch::constant_clone(CoeffType(),0.25) );
            _denom[r][s] = ant_sqrt(8.0L * (1.0L + Ms / Mr));
          }
      }

    return;
  }

  template<typename Diffusion, typename Viscosity, typename ThermalConduction, typename Mixture, typename ThermoEvaluator, typename CoeffType>
  inline
  WilkeMixture<Diffusion,Viscosity,ThermalConduction, Mixture, ThermoEvaluator,CoeffType>::~WilkeMixture()
  {
    return;
  }

  
  template<typename Diffusion, typename Viscosity, typename ThermalConduction, typename Mixture, typename ThermoEvaluator, typename CoeffType>
  inline
  CoeffType WilkeMixture<Diffusion,Viscosity,ThermalConduction, Mixture, ThermoEvaluator,CoeffType>::Mr_Ms_to_the_one_fourth( const unsigned int r,
                                                              const unsigned int s ) const
  {
    return _Mr_Ms_to_the_one_fourth[r][s];
  }
    
  
  template<typename Diffusion, typename Viscosity, typename ThermalConduction, typename Mixture, typename ThermoEvaluator, typename CoeffType>
  inline
  CoeffType WilkeMixture<Diffusion,Viscosity,ThermalConduction, Mixture, ThermoEvaluator,CoeffType>::denominator( const unsigned int r,
                                                  const unsigned int s ) const
  {
    return _denom[r][s];
  }

  template<typename Diffusion, typename Viscosity, typename ThermalConduction, typename Mixture, typename ThermoEvaluator, typename CoeffType>
  inline
  const Mixture & WilkeMixture<Diffusion,Viscosity,ThermalConduction, Mixture, ThermoEvaluator,CoeffType>::transport_mixture() const
  {
    return _transport_mixture;
  }

  template<typename Diffusion, typename Viscosity, typename ThermalConduction, typename Mixture, typename ThermoEvaluator, typename CoeffType>
  inline
  const Diffusion & WilkeMixture<Diffusion,Viscosity,ThermalConduction, Mixture, ThermoEvaluator,CoeffType>::diffusion() const
  {
     return _diffusion_set;
  }

  template<typename Diffusion, typename Viscosity, typename ThermalConduction, typename Mixture, typename ThermoEvaluator, typename CoeffType>
  inline
  const Viscosity & WilkeMixture<Diffusion,Viscosity,ThermalConduction, Mixture, ThermoEvaluator,CoeffType>::viscosity() const
  {
     return _viscosity_set;
  }

  template<typename Diffusion, typename Viscosity, typename ThermalConduction, typename Mixture, typename ThermoEvaluator, typename CoeffType>
  inline
  const ThermalConduction & WilkeMixture<Diffusion,Viscosity,ThermalConduction, Mixture, ThermoEvaluator,CoeffType>::thermal_conduction() const
  {
     return _thermal_conduction_set;
  }

  template<typename Diffusion, typename Viscosity, typename ThermalConduction, typename Mixture, typename ThermoEvaluator, typename CoeffType>
  template <typename StateType>
  inline
  const ANTIOCH_AUTO(StateType)
    WilkeMixture<Diffusion,Viscosity,ThermalConduction, Mixture, ThermoEvaluator,CoeffType>::mu(unsigned int s, const StateType & T) const
  {
     StateType mu = zero_clone(T);
     _viscosity_set(s,T,mu);
     return mu;
  }

  template<typename Diffusion, typename Viscosity, typename ThermalConduction, typename Mixture, typename ThermoEvaluator, typename CoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  void  WilkeMixture<Diffusion,Viscosity,ThermalConduction, Mixture,ThermoEvaluator,CoeffType>::mu(const StateType & T, VectorStateType & mu) const
  {
      antioch_assert_equal_to(mu.size(),_transport_mixture.n_species());
      for(unsigned int s = 0; s < _transport_mixture.n_species(); s++)
      {
          mu[s] = this->mu(s,T);
      }
  }

  template<typename Diffusion, typename Viscosity, typename ThermalConduction, typename Mixture, typename ThermoEvaluator, typename CoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  void WilkeMixture<Diffusion, Viscosity, ThermalConduction, Mixture, ThermoEvaluator, CoeffType>::D(const StateType & T, const StateType & rho,
                                                                           const StateType & cTot, const VectorStateType & mass_fractions, 
                                                                           const VectorStateType & mu, VectorStateType & ds) const
  {
// first way
// \todo: if not needed, find a way to supress Ds building
// \todo: initialization as full squared matrix, even though half is needed and used
     typename rebind<VectorStateType,VectorStateType>::type Ds(ds.size());
     init_constant(Ds,ds);
     _diffusion_set(T,cTot,Ds );
     wilke_diffusion_rule(_transport_mixture.chemical_mixture(), mass_fractions, Ds, ds, typename physical_tag<typename Diffusion::model>::diffusion_species_type ());

// second way
// \todo: if not needed, find a way to supress k building
// \todo: Ds[s][s] not needed here
     StateType k = zero_clone(T);
     for(unsigned int s = 0; s < _transport_mixture.n_species(); s++)
     {
        _thermal_conduction_set(s,mu[s],Ds[s][s],T,cTot * _transport_mixture.chemical_mixture().M(s),k);
        _diffusion_set(_thermo_evaluator.cp(TempCache<StateType>(T),s), k, ds[s]);
     }
  }

  template<typename Diffusion, typename Viscosity, typename ThermalConduction, typename Mixture, typename ThermoEvaluator, typename CoeffType>
  template <typename StateType>
  inline
  const ANTIOCH_AUTO(StateType) 
        WilkeMixture<Diffusion,Viscosity,ThermalConduction, Mixture, ThermoEvaluator,CoeffType>::k(unsigned int s, const StateType & mu, 
                                                                                                   const StateType & T, const StateType &n_molar_mixture) const
  {
// if needed
// \todo: if not needed, find a way to supress Ds building
     StateType Ds = zero_clone(T);
     _diffusion_set(s,T,n_molar_mixture,Ds);

     StateType k = zero_clone(T);
     _thermal_conduction_set(s,mu,Ds,T,n_molar_mixture * _transport_mixture.chemical_mixture().M(s),k);

     return k;
  }

  template<typename Diffusion, typename Viscosity, typename ThermalConduction, typename Mixture, typename ThermoEvaluator, typename CoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  void WilkeMixture<Diffusion,Viscosity,ThermalConduction, Mixture, ThermoEvaluator,CoeffType>::k(const VectorStateType & mu, const StateType & T, 
                                                                                                  const VectorStateType & mass_fractions, 
                                                                                                  const StateType &rho, VectorStateType & k) const
  {
      antioch_assert_equal_to(k.size(),mu.size());
      antioch_assert_equal_to(k.size(),mass_fractions.size());
      antioch_assert_equal_to(k.size(),_transport_mixture.n_species());

    const StateType n_molar_mixture = rho / _transport_mixture.chemical_mixture().M(mass_fractions);
     for(unsigned int s = 0; s < _transport_mixture.n_species(); s++)
     {
          k[s] = this->k(s,mu[s],T,n_molar_mixture);
     }
  }


  // species mu and k
  // mixture D
  template<typename Diffusion, typename Viscosity, typename ThermalConduction, typename Mixture, typename ThermoEvaluator, typename CoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  void WilkeMixture<Diffusion, Viscosity, ThermalConduction, Mixture, ThermoEvaluator, CoeffType>::D_and_k(const VectorStateType & mu, const StateType & T, const StateType & rho, 
                                                                                                           const VectorStateType & mass_fractions, 
                                                                                                           VectorStateType & k, VectorStateType & ds) const
  {

     antioch_assert_equal_to(ds.size(),mass_fractions.size());
     antioch_assert_equal_to(ds.size(),_transport_mixture.n_species());
     antioch_assert_equal_to(mu.size(),_transport_mixture.n_species());
     const StateType n_molar_mixture = rho / _transport_mixture.chemical_mixture().M(mass_fractions); // total molar density

// diffusion comes first
// \todo: if not needed, find a way to supress Ds building
// \todo: initialization as full squared matrix, even though half is needed and used
     typename Antioch::rebind<VectorStateType,VectorStateType>::type Ds(mass_fractions.size());
     init_constant(Ds,ds);
     _diffusion_set(T, n_molar_mixture, Ds);
     wilke_diffusion_rule(_transport_mixture.chemical_mixture(), mass_fractions, Ds, ds, typename physical_tag<typename Diffusion::model>::diffusion_species_type ());

// thermal conduction
    for(unsigned int s = 0; s < _transport_mixture.n_species(); s++)
    {
      _diffusion_set(s,T,n_molar_mixture,Ds[s][s]);
      _thermal_conduction_set(s,mu[s],Ds[s][s],T,n_molar_mixture * _transport_mixture.chemical_mixture().M(s),k[s]);
    }

// diffusion comes last
    for(unsigned int s = 0; s < _transport_mixture.n_species(); s++)
        _diffusion_set(rho, _thermo_evaluator.cp(TempCache<StateType>(T),s), k[s], ds[s] );
  }


} // end namespace Antioch

#endif // ANTIOCH_WILKE_MIXTURE_H
