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

// C++
#include <vector>

namespace Antioch
{
// Forward declarations

  template<class CoeffType>
  class TransportMixture;

  template<class CoeffType>
  class TempCache;

  template<class Diffusion, class Viscosity, class ThermalConduction,
           class CoeffType=double>
  class WilkeMixture
  {
  public:

    WilkeMixture( const Diffusion & diffusion, const Viscosity & viscosity, 
                  const ThermalConduction & thermal_conduction,
                  const TransportMixture<CoeffType>& transport_mixture );
    ~WilkeMixture();

    const TransportMixture<CoeffType>& transport_mixture() const; // contains the thermo

    const PhysicalSet<Diffusion> & diffusion() const;

    const PhysicalSet<Viscosity> & viscosity() const;

    const PhysicalSet<ThermalConduction> & thermal_conduction() const;

    template <typename StateType, typename VectorStateType>
    void D(unsigned int s, const StateType & T, const StateType & cTot, const VectorStateType & mass_fractions, const StateType mu, VectorStateType & ds) const;

    template <typename StateType>
    const ANTIOCH_AUTO(StateType)
     mu(unsigned int s, const StateType & T) const
     ANTIOCH_AUTOFUNC(StateType, _viscosity_set(s,T,typename physical_tag<Viscosity>::type))

    template <typename StateType, typename VectorStateType>
    void  mu(const StateType & T, VectorStateType & mu) const;

    template <typename StateType, typename VectorStateType>
    const StateType k(unsigned int s, const StateType & mu, const StateType & T, const StateType &rho, const VectorStateType & mass_fractions) const;

    template <typename StateType, Typename VectorStateType>
    void D_and_k(const StateType & mu, unsigned int s, const StateType & T, const VectorStateType & mass_fractions, 
                       VectorStateType & k, VectorStateType & D) const;

  protected:

    CoeffType Mr_Ms_to_the_one_fourth( const unsigned int r,
                                       const unsigned int s ) const;
    
    CoeffType denominator( const unsigned int r,
                           const unsigned int s ) const;

    const TransportMixture<CoeffType>& _transport_mixture;

    const PhysicalSet<Diffusion>         & _diffusion_set;

    const PhysicalSet<Viscosity>         & _viscosity_set;

    const PhysicalSet<ThermalConduction> & _thermal_conduction_set;

    //! Cache for numerator term
    /*! \todo We should use a more efficient data structure */
    std::vector<std::vector<CoeffType> > _Mr_Ms_to_the_one_fourth;
    
    //! Cache for denominator term
    /*! \todo We should use a more efficient data structure */
    std::vector<std::vector<CoeffType> > _denom;

  };

  template<class Diffusion, class Viscosity, class ThermalConduction, class CoeffType>
  WilkeMixture<Diffusion,Viscosity,ThermalConduction,CoeffType>::WilkeMixture(const Diffusion & diffusion, const Viscosity & viscosity, 
                                        const ThermalConduction & thermal_conduction,
                                        const TransportMixture<CoeffType>& transport_mixture )
    : _transport_mixture(transport_mixture),
      _diffusion_set(diffusion),
      _viscosity_set(viscosity),
      _thermal_conduction_set(thermal_conduction),
      _Mr_Ms_to_the_one_fourth(_transport_mixture.n_species(), std::vector<CoeffType>(_transport_mixture.n_species())),
      _denom(_transport_mixture.n_species(), std::vector<CoeffType>(_transport_mixture.n_species()))
  {

    for( unsigned int r = 0; r < transport_mixture.n_species(); r++ )
      {
        for( unsigned int s = 0; s < transport_mixture.n_species(); s++ )
          {
            const CoeffType Mr = _transport_mixture.M(r);
            const CoeffType Ms = _transport_mixture.M(s);

            _Mr_Ms_to_the_one_fourth[r][s] = ant_pow( Mr/Ms, Antioch::constant_clone(CoeffType,0.25) );
            _denom[r][s] = ant_sqrt(8.0L * (1.0L + Ms / Mr));
          }
      }

    return;
  }

  template<class Diffusion, class Viscosity, class ThermalConduction, class CoeffType>
  WilkeMixture<Diffusion,Viscosity,ThermalConduction,CoeffType>::~WilkeMixture()
  {
    return;
  }

  
  template<class Diffusion, class Viscosity, class ThermalConduction, class CoeffType>
  inline
  CoeffType WilkeMixture<Diffusion,Viscosity,ThermalConduction,CoeffType>::Mr_Ms_to_the_one_fourth( const unsigned int r,
                                                              const unsigned int s ) const
  {
    return _Mr_Ms_to_the_one_fourth[r][s];
  }
    
  
  template<class Diffusion, class Viscosity, class ThermalConduction, class CoeffType>
  inline
  CoeffType WilkeMixture<Diffusion,Viscosity,ThermalConduction,CoeffType>::denominator( const unsigned int r,
                                                  const unsigned int s ) const
  {
    return _denom[r][s];
  }

  template<class Diffusion, class Viscosity, class ThermalConduction, class CoeffType>
  template <typename StateType>
  inline
  const ANTIOCH_AUTO(StateType) WilkeMixture<Diffusion,Viscosity,ThermalConduction,CoeffType>::k_trans(unsigned int s, const StateType & T, const StateType &rho) const
  {
     
  }

  template<class Diffusion, class Viscosity, class ThermalConduction, class CoeffType>
  inline
  const TransportMixture<CoeffType>& WilkeMixture<Diffusion,Viscosity,ThermalConduction,CoeffType>::transport_mixture() const
  {
    return _transport_mixture;
  }

  template<class Diffusion, class Viscosity, class ThermalConduction, class CoeffType>
  inline
  const PhysicalSet<Diffusion> & WilkeMixture<Diffusion,Viscosity,ThermalConduction,CoeffType>::diffusion() const
  {
     return _diffusion;
  }

  template<class Diffusion, class Viscosity, class ThermalConduction, class CoeffType>
  inline
  const PhysicalSet<Viscosity> & WilkeMixture<Diffusion,Viscosity,ThermalConduction,CoeffType>::viscosity() const
  {
     return _viscosity;
  }

  template<class Diffusion, class Viscosity, class ThermalConduction, class CoeffType>
  inline
  const PhysicalSet<ThermalConduction> & WilkeMixture<Diffusion,Viscosity,ThermalConduction,CoeffType>::thermal_conduction() const
  {
     return _thermal_conduction;
  }

  template<class Diffusion, class Viscosity, class ThermalConduction, class CoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  void  WilkeMixture<Diffusion,Viscosity,ThermalConduction,CoeffType>::mu(const StateType & T, VectorStateType & mu) const
  {
      antioch_assert_equal_to(mu.size(),_transport_mixture.n_species());
      for(unsigned int s = 0; s < _transport_mixture.n_species(); s++)
      {
          mu[s] = this->mu(s,T);
      }
  }

  template<class Diffusion, class Viscosity, class ThermalConduction, class CoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  void WilkeMixture<Diffusion, Viscosity, ThermalConduction, CoeffType>::D(const StateType & T, const StateType & cTot, const VectorStateType & mass_fractions, 
                                                                           const StateType & mu, VectorStateType & ds) const
  {
// first way
// TODO: if not needed, find a way to supress Ds building
     typename rebind<VectorStateType,VectorStateType>::type Ds(ds.size(),zero_clone(ds));
     _diffusion_set(T,cTot,Ds, typename physical_tag<Diffusion>::type());
     wilke_diffusion_rule(_transport_mixture.chemical_mixture(),mass_fractions,Ds,ds,physical_tag<Diffusion>::type());
// second way
// TODO: if not needed, find a way to supress k building
     StateType k = zero_clone(cTot);
     for(unsigned int s = 0; s < _transport_mixture.n_species(); s++)
     {
        _thermal_conduction_set(s,mu,D[s][s],T,rho,k, typename physical_tag<ThermalConduction>::type());
        _diffusion_set(_transport_mixture.thermo().cp(TempCache(T),s), k[s], ds[s], typename physical_tag<Diffusion>::type());
     }
  }

  template<class Diffusion, class Viscosity, class ThermalConduction, class CoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  const StateType WilkeMixture<Diffusion,Viscosity,ThermalConduction,CoeffType>::k(unsigned int s, const StateType & mu, const StateType & T, const StateType &rho, const VectorStateType & mass_fractions) const
  {
     StateType k = zero_clone(T);
// if needed
// TODO: if not needed, find a way to supress Ds building
     VectorStateType Ds = zero_clone(mass_fractions);
     _diffusion_set(T,rho / _transport_mixture().chemical_mixture().M(mass_fractions),Ds, typename physical_tag<Diffusion>::second_type());

     _thermal_conduction_set(s,mu,D[s][s],T,rho,k, typename physical_tag<ThermalConduction>::type());

     return k;
  }


  // species mu and k
  // mixture D
  template <typename StateType, typename VectorStateType>
  inline
  void WilkeMixture<Diffusion, Viscosity, ThermalConduction, CoeffType>::D_and_k(const VectorStateType & mu, const StateType & T, const StateType & rho, const VectorStateType & mass_fractions, 
                                                                                       VectorStateType & k, VectorStateType & ds) const
  {

     antioch_assert_equal_to(ds.size(),mass_fractions.size());
     antioch_assert_equal_to(ds.size(),_transport_mixture.n_species());
     antioch_assert_equal_to(mu.size(),_transport_mixture.n_species());

// diffusion comes first
// TODO: if not needed, find a way to supress Ds building and cTot computation
     VectorStateType Ds = zero_clone(mass_fractions);
     _diffusion_set(T, rho / _transport_mixture().chemical_mixture().M(mass_fractions), Ds, typename physical_tag<Diffusion>::type());
     wilke_diffusion_rule(_transport_mixture.chemical_mixture(), mass_fractions, Ds, ds, physical_tag<Diffusion>::type());

// thermal conduction
    for(unsigned int s = 0; s < _transport_mixture.n_species(); s++)
        _thermal_conduction_set(s,mu[s],Ds[s][s],T,rho,k[s], typename physical_tag<ThermalConduction>::type());

// diffusion comes last
    for(unsigned int s = 0; s < _transport_mixture.n_species(); s++)
        _diffusion_set(_transport_mixture.thermo().cp(TempCache(T),s), k[s], ds[s], typename physical_tag<Diffusion>::type());
  }


} // end namespace Antioch

#endif // ANTIOCH_WILKE_MIXTURE_H
