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
//--------------------------------------------------------------------------

#ifndef ANTIOCH_KINETICS_THEORY_VISCOSITY_H
#define ANTIOCH_KINETICS_THEORY_VISCOSITY_H

// Antioch
#include "antioch/metaprogramming.h"
#include "antioch/antioch_asserts.h"
#include "antioch/viscosity_enum.h"
#include "antioch/math_constants.h"
#include "antioch/cmath_shims.h"
#include "antioch/chemical_mixture.h"
#include "antioch/Stockmayer_potential.h"
#include "antioch/Lennard_Jones_potential.h"

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
  template<typename CoeffType = double >
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
       * in the calculation of _delta_star,
       * we use the coefficient to SI from debye (3.335641 10-30 from http://cccbdb.nist.gov/debye.asp)
       */
      PureSpeciesViscosity(const CoeffType & LJ_depth, const CoeffType & LJ_diameter, // depth in K (epsilon/kB), diameter in angström
                               const CoeffType & dipole_moment, const CoeffType & mass);  // dipole moment in D, molecular mass in kg 

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
      ANTIOCH_AUTOFUNC(StateType,    _a * CoeffType(1e-13)  // 5 / 16 * sqrt(pi * Boltzmann_constant)
                                        * ant_sqrt(CoeffType(1e26L) * _mass * T )  
                                     / ( Constants::pi<CoeffType>() * _LJ.diameter() * _LJ.diameter() * CoeffType(1e-20L) * // to SI
                                          collision_integral( T )   // Omega(2,2), T*
                                        )
                      )

      template <typename StateType>
      ANTIOCH_AUTO(StateType) 
      derivative(const StateType &T) const
      ANTIOCH_AUTOFUNC(StateType,  this->viscosity(T) * 
                           (StateType (1.L)/T - dcollision_integral_dT(T) / collision_integral(T)     // T*
                           ))

      template <typename StateType>
      StateType compute_viscosity_and_derivative( const StateType& T, StateType & viscosity, StateType & dviscosity_dT ) const;

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

      /*! never ever use it*/
      PureSpeciesViscosity();

      //! angles integrated collision integral
      template <typename StateType>
      const StateType collision_integral(const StateType & T) const;

      //! angles integrated collision integral
      template <typename StateType>
      const StateType dcollision_integral_dTstar(const StateType & T) const;

      const CoeffType _a;

      LennardJonesPotential<CoeffType> _LJ;
      CoeffType                        _dipole_moment;
      CoeffType                        _mass;

      CoeffType              _delta_star;
      std::vector<CoeffType> _collision_coefficients;
      std::vector<CoeffType> _temperature_line;
      std::vector<CoeffType> _temperature;

      const ViscosityModel _model;

  };


  template <typename CoeffType>
  inline
  PureSpeciesViscosity<CoeffType>::PureSpeciesViscosity(const CoeffType & LJ_depth, const CoeffType & LJ_diameter, 
                                                              const CoeffType & dipole_moment, const CoeffType & mass):
        Viscosity<CoeffType>(ViscosityModel::KINETICSTHEORY),
        _a(0.3125e-12L * ant_sqrt(Constants::pi<CoeffType>() * Constants::Boltzmann_constant<CoeffType>() * 1e24L)), /* 5 / 16 * sqrt(pi * Boltzmann constant) */
        _LJ(LJ_depth,LJ_diameter),
        _dipole_moment(dipole_moment),
        _mass(mass),
        _delta_star(ant_pow(_dipole_moment * CoeffType(3.335641L),2) * CoeffType(1e-30L) /             
                     ( _LJ.depth() * CoeffType(8.L) * Constants::pi<CoeffType>() * Constants::vacuum_permittivity<CoeffType>() * 
                           Constants::Boltzmann_constant<CoeffType>() * ant_pow(_LJ.diameter(),3) )),
        _collision_coefficients(3,0.),  // from ChemKin: quadratic interpolation
        _model(PURE_SPECIES)
  {
     StockmayerPotential<CoeffType> stock;
     stock.temperature_interpolation(StockmayerPotential<CoeffType>::INTEGRAL::TWO,_delta_star,_collision_coefficients,_temperature_line);
     _temperature = stock.temperature();

     return;
  }


  template <typename CoeffType>
  inline
  PureSpeciesViscosity<CoeffType>::~PureSpeciesViscosity()
  {
     return;
  }

  template <typename CoeffType>
  inline
  void PureSpeciesViscosity<CoeffType>::reset_coeffs( const CoeffType & LJ_depth, const CoeffType & LJ_dia, const CoeffType & dipole_moment, const CoeffType & mass )
  {
//redefining parameters
     _LJ.reset_coeffs(LJ_depth,LJ_dia);
     _dipole_moment = dipole_moment;
     _mass = mass;
     _delta_star = (ant_pow(_dipole_moment * CoeffType(3.335641L),2) * CoeffType(1e-30L) /             
                     ( _LJ.depth() * CoeffType(8.L) * Constants::pi<CoeffType>() * Constants::vacuum_permittivity<CoeffType>() * 
                           Constants::Boltzmann_constant<CoeffType>() * ant_pow(_LJ.diameter(),3) ));

//redefining collision integral
     _collision_coefficients.clear();
     _collision_coefficients.resize(3,0.);

     StockmayerPotential<CoeffType> stock;
     stock.temperature_interpolation(StockmayerPotential<CoeffType>::INTEGRAL::TWO,_delta_star,_collision_coefficients,_temperature_line);
  }

  template <typename CoeffType>
  inline
  void PureSpeciesViscosity<CoeffType>::reset_coeffs( const std::vector<CoeffType> & coeffs )
  {
     antioch_assert_equal_to(coeffs.size(),4);

     this->reset_coeffs(coeffs[0],coeffs[1],coeffs[2],coeffs[3]);
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  const StateType PureSpeciesViscosity<CoeffType>::collision_integral(const StateType & T) const
  {
     StateType omega_22 = zero_clone(T);
     StateType Tmp = constant_clone(T,1);

     // Stockmayer potential integration on log of temperature
     StateType Tlogstar = ant_log(T / _LJ.depth() );
     std::vector<StateType> coeffs(_collision_coefficients.size());

     PolynomialRegression reg;

     StateType eps =  reg.polynomial_regression(_temperature,_temperature_line,coeffs,_collision_coefficients.size(),Tlogstar);
     if(min(eps) < 0)
     {
        std::cerr << "Failed regression: " << eps << std::endl;
        antioch_error();
     }
     

     for(unsigned int d = 0; d < _collision_coefficients.size(); d++)
     {
        omega_22 += Tmp * coeffs[d];
        Tmp = Tmp * Tlogstar;
     }
     return omega_22;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  const StateType PureSpeciesViscosity<CoeffType>::dcollision_integral_dTstar(const StateType & T) const
  {
     StateType domega22_dTstar = zero_clone(T);
     StateType Tmp = constant_clone(T,1);
     StateType Tlogstar = ant_log(T / _LJ.depth() );

     for(unsigned int d = 1; d < _collision_coefficients.size(); d++)
     {
        domega22_dTstar += constant_clone(T,d) * Tmp * _collision_coefficients[d];
        Tmp *= Tlogstar;
     }
     return domega22_dTstar;
  }



  template <typename CoeffType>
  template <typename StateType>
  inline
  StateType PureSpeciesViscosity<CoeffType>::compute_viscosity_and_derivative( const StateType& T, StateType & viscosity, StateType & dviscosity_dT ) const
  {
     viscosity = this->viscosity(T);
     dviscosity_dT = viscosity *
                           (StateType (1.)/T - dcollision_integrale_dT( T / _LJ.depth() ) /(_LJ.depth() * collision_integrale( T / _LJ.depth()) ));  // T*, dc/dT = dc/dT* * dT*/dT = dc/dT* / _LJ.depth() 
     return;
  }

  template <typename CoeffType>
  inline
  void PureSpeciesViscosity<CoeffType>::print(std::ostream& os) const
  {
     os << "Kinetics theory viscosity:\n"
        << "5/16 * sqrt(pi * " << _mass << " * kb * T) / ( pi * " << _LJ.depth() << "^2 * <Omega(2,2)*> )\n"
        << "T* = T / " << _LJ.depth() << "\n"
        << "delta* = 1/2 * " << _dipole_moment << "^2 / ( " << _LJ.depth() << " * kb * " << _LJ.diameter() << "^3 )";
  }

  template <typename CoeffType>
  inline
  const CoeffType & PureSpeciesViscosity<CoeffType>::delta_star() const
  {
     return _delta_star;
  }

  template<typename CoeffType>
  inline
  ViscosityModel PureSpeciesViscosity<CoeffType>::model() const
  {
      return _model;
  }

} // end namespace Antioch

#endif //ANTIOCH_BLOTTNER_VISCOSITY_H
