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
//--------------------------------------------------------------------------

#ifndef ANTIOCH_PURE_SPECIES_VISCOSITY_H
#define ANTIOCH_PURE_SPECIES_VISCOSITY_H

// Antioch
#include "antioch/metaprogramming.h"
#include "antioch/antioch_asserts.h"
#include "antioch/math_constants.h"
#include "antioch/cmath_shims.h"
#include "antioch/chemical_mixture.h"
#include "antioch/Stockmayer_potential.h"
#include "antioch/gsl_spliner.h"
#include "antioch/Lennard_Jones_potential.h"
#include "antioch/physical_constants.h"

// C++
#include <cmath>
#include <vector>
#include <iostream>

namespace Antioch
{

  /*! \class PureSpeciesViscosity
   *
   * Get the kinetics theory's expression of
   * viscosity. Ideal gaz assumption.
   *
   * \f[
   *     \nu_k = \frac{5}{16} \frac{\sqrt{\pi \mathrm{k_B} T}}{\pi \sigma_k^2 \Omega^{(2,2)*}}
   * \f]
   * with:
   *  - \f$\mathrm{k_B}\f$ the Boltzmann constant
   *  - \f$\sigma_k\f$ the Lennard-Jones collision diameter of molecule \f$k\f$
   *  - \f$\Omega^{(2,2)*}\f$ the collision integral, determined by a fit using: 
   *            - \f[
   *                  T^* = \frac{\mathrm{k_B} T}{\epsilon_k}
   *              \f]
   *            - \f[
   *                  \delta_k^* = \frac{1}{2} \frac{\alpha_k^2}{\epsilon_k\sigma_k^3}
   *              \f]
   * with \f$\epsilon\f$ the Lennard-Jones potential well depth and \f$\alpha\f$ the
   * dipole moment.
   */
  template<typename CoeffType = double, typename Interpolator = GSLSpliner>
  class PureSpeciesViscosity
  {
    public:

      /*!\brief Constructor
       *
       * It takes:
       *   - Lennard-Jones depth in K (depth over Boltzmann constant)
       *   - Lennard-Jones diameter in angström
       *   - dipole moment in debye
       *   - mass in kg
       *
       * Angström,  debye and Boltzmann constant will (almost) cancel out the orders of magnitude
       * in the calculation of _delta_star.
       */
      PureSpeciesViscosity(const CoeffType & LJ_depth, const CoeffType & LJ_diameter, // depth in K (epsilon/kB), diameter in angström
                               const CoeffType & dipole_moment, const CoeffType & mass);  // dipole moment in D, molecular mass in kg 
      PureSpeciesViscosity(const std::vector<CoeffType> & coeffs);

      ~PureSpeciesViscosity();


      void reset_coeffs( const CoeffType & LJ_depth, const CoeffType & LJ_dia, const CoeffType & dipole_moment, const CoeffType & mass );
      void reset_coeffs( const std::vector<CoeffType> & coeffs );

      template <typename StateType>
      ANTIOCH_AUTO(StateType) 
      operator()(const StateType &T) const
      ANTIOCH_AUTOFUNC(StateType,  this->viscosity(T)  )

      /*! \brief Calculates kinetics theory viscosity
       *
       * input is temperature (K), uses
       *   - Lennard-Jones diameter (in Angström)
       *   - reduced Lennard-Jones depth (in K)
       *   - mass of the molecule (in kg)
       *   - integrated collision integral given by the
       *        Stockmayer potential (no unit)
       */
      template <typename StateType>
      ANTIOCH_AUTO(StateType) 
      viscosity(const StateType &T) const
      ANTIOCH_AUTOFUNC(StateType,    _a   // 5 / 16 * sqrt(pi * Boltzmann_constant)
                                        * ant_sqrt(_mass * T )  
                                     / ( Constants::pi<CoeffType>() * _LJ.diameter() * _LJ.diameter() * CoeffType(1e-20L) * // to SI
                                         _interp.interpolated_value( T / _LJ.depth() )   // Omega(2,2), T*
                                        )
                      )

      template <typename StateType>
      ANTIOCH_AUTO(StateType) 
      derivative(const StateType &T) const
      ANTIOCH_AUTOFUNC(StateType,  this->viscosity(T) * 
                           (StateType (1.L)/T - _interp.dinterp_dx(T / _LJ.depth()) / _interp.interpolated_value(T / _LJ.depth())     // T*
                           ))

      template <typename StateType>
      StateType compute_viscosity_and_derivative( const StateType& T, StateType & viscosity, StateType & dviscosity_dT ) const;

      template <typename StateType>
      ANTIOCH_AUTO(StateType)
        Stockmayer(const StateType & T) const
      ANTIOCH_AUTOFUNC(StateType,_interp.interpolated_value( T / _LJ.depth() ) )   // Omega(2,2)
                                          

      //! Formatted print, by default to \p std::cout
      void print(std::ostream& os = std::cout) const;

      //!\return the value of the reduced dipole moment
      const CoeffType & delta_star() const;

      //! Formatted print.
      friend std::ostream& operator<<(std::ostream& os, const PureSpeciesViscosity& mu)
      {
        mu.print(os);
        return os;
      }
      

    private:

      void build_interpolation();

      /*! never ever use it*/
      PureSpeciesViscosity();

      const CoeffType _a;

      LennardJonesPotential<CoeffType> _LJ;
      CoeffType                        _dipole_moment;
      CoeffType                        _mass;


      Interpolator    _interp;
      CoeffType _delta_star;

  };

  template <typename CoeffType, typename Interpolator>
  inline
  PureSpeciesViscosity<CoeffType,Interpolator>::PureSpeciesViscosity(const CoeffType & LJ_depth, const CoeffType & LJ_diameter, 
                                                              const CoeffType & dipole_moment, const CoeffType & mass):
        _a(0.3125L * ant_sqrt(Constants::pi<CoeffType>() * Constants::Boltzmann_constant<CoeffType>())), /* 5 / 16 * sqrt(pi * Boltzmann constant) */
        _LJ(LJ_depth,LJ_diameter),
        _dipole_moment(dipole_moment),
        _mass(mass),
        _delta_star(ant_pow(_dipole_moment * Units<CoeffType>("D").get_SI_factor(),2) /             
                     ( _LJ.depth() * CoeffType(8.L) * Constants::pi<CoeffType>() * Constants::vacuum_permittivity<CoeffType>() * 
                           Constants::Boltzmann_constant<CoeffType>() * ant_pow(_LJ.diameter(),3) ))
  {
std::cout << std::setprecision(15) << Constants::pi<CoeffType>() << " * " << Constants::Boltzmann_constant<CoeffType>() << std::endl;
     this->build_interpolation();
     return;
  }

  template <typename CoeffType, typename Interpolator>
  inline
  PureSpeciesViscosity<CoeffType,Interpolator>::PureSpeciesViscosity(const std::vector<CoeffType> & coeffs):
        _a(0.3125 * ant_sqrt(Constants::pi<CoeffType>() * Constants::Boltzmann_constant<CoeffType>())), /* 5 / 16 * sqrt(pi * Boltzmann constant) */
#ifndef NDEBUG
        _LJ(-1,-1),
        _dipole_moment(-1),
        _mass(-1),
        _delta_star(-1)
#else
        _LJ(coeffs[0],coeffs[1]),
        _dipole_moment(coeffs[2]),
        _mass(coeffs[3]),
        _delta_star(ant_pow(_dipole_moment * Units<CoeffType>("D").get_SI_factor(),2) /
                     ( _LJ.depth() * CoeffType(8.L) * Constants::pi<CoeffType>() * Constants::vacuum_permittivity<CoeffType>() * 
                           Constants::Boltzmann_constant<CoeffType>() * ant_pow(_LJ.diameter(),3) ))
#endif
  {
#ifndef NDEBUG
        antioch_assert_equal_to(coeffs.size(),4);

        _LJ.set_depth(coeffs[0]);
        _LJ.set_diameter(coeffs[1]);
        _dipole_moment(coeffs[2]);
        _mass(coeffs[3]);
        _delta_star = ant_pow(_dipole_moment * Units<CoeffType>("D").get_SI_factor(),2) /
                     ( _LJ.depth() * CoeffType(8.L) * Constants::pi<CoeffType>() * Constants::vacuum_permittivity<CoeffType>() * 
                           Constants::Boltzmann_constant<CoeffType>() * ant_pow(_LJ.diameter(),3) );
#endif
     this->build_interpolation();
     return;
  }


  template <typename CoeffType, typename Interpolator>
  inline
  PureSpeciesViscosity<CoeffType,Interpolator>::~PureSpeciesViscosity()
  {
     return;
  }

  template <typename CoeffType, typename Interpolator>
  inline
  void PureSpeciesViscosity<CoeffType,Interpolator>::build_interpolation()
  {

     _interp.spline_delete();
     StockmayerPotential<CoeffType> surface;
     std::vector<CoeffType> interp_surf(surface.temperature().size(),0);
     for(unsigned int iT = 0; iT < surface.temperature().size(); iT++)
     {
        Interpolator spline(surface.delta(),surface.omega_2_2()[iT]);
        interp_surf[iT] = spline.interpolated_value(_delta_star);
     }

     _interp.spline_init(surface.temperature(),interp_surf);
  }

  template <typename CoeffType, typename Interpolator>
  inline
  void PureSpeciesViscosity<CoeffType,Interpolator>::reset_coeffs( const CoeffType & LJ_depth, const CoeffType & LJ_dia, const CoeffType & dipole_moment, const CoeffType & mass )
  {
//redefining parameters
     _LJ.reset_coeffs(LJ_depth,LJ_dia);
     _dipole_moment = dipole_moment;
     _mass = mass;
     _delta_star = (ant_pow(_dipole_moment * Units<CoeffType>("D").get_SI_factor(),2) /
                     ( _LJ.depth() * CoeffType(8.L) * Constants::pi<CoeffType>() * Constants::vacuum_permittivity<CoeffType>() * 
                           Constants::Boltzmann_constant<CoeffType>() * ant_pow(_LJ.diameter(),3) ));

//redefining collision integral
    this->build_interpolation();
  }

  template <typename CoeffType, typename Interpolator>
  inline
  void PureSpeciesViscosity<CoeffType,Interpolator>::reset_coeffs( const std::vector<CoeffType> & coeffs )
  {
     antioch_assert_equal_to(coeffs.size(),4);

     this->reset_coeffs(coeffs[0],coeffs[1],coeffs[2],coeffs[3]);
  }

  template <typename CoeffType, typename Interpolator>
  template <typename StateType>
  inline
  StateType PureSpeciesViscosity<CoeffType,Interpolator>::compute_viscosity_and_derivative( const StateType& T, StateType & viscosity, StateType & dviscosity_dT ) const
  {
     viscosity = this->viscosity(T);
     dviscosity_dT = viscosity *
                           (StateType (1.)/T - _interp.dinterp_dx( T / _LJ.depth() ) /(_LJ.depth() * _interp.interpolated_value( T / _LJ.depth()) ));  // T*, dc/dT = dc/dT* * dT*/dT = dc/dT* / _LJ.depth() 
     return;
  }

  template <typename CoeffType, typename Interpolator>
  inline
  void PureSpeciesViscosity<CoeffType,Interpolator>::print(std::ostream& os) const
  {
     os << "Pure species viscosity:\n"
        << "5/16 * sqrt(pi * " << _mass << " * kb * T) / ( pi * " << _LJ.depth() << "^2 * <Omega(2,2)*> )\n"
        << "T* = T / " << _LJ.depth() << "\n"
        << "delta* = 1/2 * " << _dipole_moment << "^2 / ( " << _LJ.depth() << " * kb * " << _LJ.diameter() << "^3 )";
  }

  template <typename CoeffType, typename Interpolator>
  inline
  const CoeffType & PureSpeciesViscosity<CoeffType,Interpolator>::delta_star() const
  {
     return _delta_star;
  }

} // end namespace Antioch

#include "antioch/pure_species_viscosity_utils_decl.h"

#endif //ANTIOCH_PURE_SPECIES_VISCOSITY_H
