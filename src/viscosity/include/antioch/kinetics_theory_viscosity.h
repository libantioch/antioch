//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014-2016 Paul T. Bauman, Benjamin S. Kirk,
//                         Sylvain Plessis, Roy H. Stonger
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

#include "antioch_config.h"
#ifdef ANTIOCH_HAVE_GSL // if we do not have it, we don't even define the stuff

#ifndef ANTIOCH_KINETICS_THEORY_VISCOSITY_H
#define ANTIOCH_KINETICS_THEORY_VISCOSITY_H

// Antioch
#include "antioch/metaprogramming.h"
#include "antioch/antioch_asserts.h"
#include "antioch/math_constants.h"
#include "antioch/cmath_shims.h"
#include "antioch/chemical_mixture.h"
#include "antioch/stockmayer_potential.h"
#include "antioch/gsl_spliner.h"
#include "antioch/lennard_jones_potential.h"
#include "antioch/physical_constants.h"
#include "antioch/species_viscosity_base.h"

// C++
#include <cmath>
#include <vector>
#include <iostream>

namespace Antioch
{

  /*! \class KineticsTheoryViscosity
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
   *
   * An important remark about precision:
   * the viscosity is computed as \f$c * s(T)\f$ with \f$c\f$ a coefficient
   * and \f$s(T)\f$ the Stockmayer potential spline. The expression of \f$c\f$
   * is
   * \f[
   *     c = \frac{5}{16} \sqrt{\frac{\mathrm{k_B} m}{\pi}} \frac{1}{\sigma^2}
   * \f]
   *  which roughly scales as
   * \f[
   *    10^{-6}  = 10^{-1} \sqrt{\frac{10^{-23} 10^{-27}}{1}} \frac{1}{10^{-20}}
   * \f]
   * The main issue lies in the square root, the float precision evaluate
   * \f$\sqrt{10^{-50}}\f$ to be zero. Therefore some adaptation has been
   * made by multiplying the factor in the square root by \f$10^{28}\f$ and
   * multiplying afterwards by \f$10^{-14}\f$.
   */
  template<typename CoeffType = double, typename Interpolator = GSLSpliner>
  class KineticsTheoryViscosity : public SpeciesViscosityBase<KineticsTheoryViscosity<CoeffType,Interpolator>,CoeffType>
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

      KineticsTheoryViscosity(const CoeffType & LJ_depth, // depth in K (epsilon/kB),
                              const CoeffType & LJ_diameter, // diameter in angström
                              const CoeffType & dipole_moment,
                              const CoeffType & mass);

      // dipole moment in D, molecular mass in kg
      KineticsTheoryViscosity(const std::vector<CoeffType> & coeffs);

      virtual ~KineticsTheoryViscosity(){};


      void reset_coeffs( const CoeffType & LJ_depth, const CoeffType & LJ_dia, const CoeffType & dipole_moment, const CoeffType & mass );



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
      ANTIOCH_AUTOFUNC(StateType,    _a   // 5 / 16 * sqrt(pi * Boltzmann_constant * mass) / ( pi * sigma * sigma )
                                       * _interp.interpolated_value(T)    // sqrt(T) / Omega<(2,2)>(T*)
                      )

      template <typename StateType>
      ANTIOCH_AUTO(StateType)
      derivative(const StateType &T) const
      ANTIOCH_AUTOFUNC(StateType,  _a * _interp.dinterp_dx(T) )

      template <typename StateType>
      StateType compute_viscosity_and_derivative( const StateType& T, StateType & viscosity, StateType & dviscosity_dT ) const;

      template <typename StateType>
      ANTIOCH_AUTO(StateType)
        Stockmayer(const StateType & T) const
      ANTIOCH_AUTOFUNC(StateType,ant_sqrt(T) / _interp.interpolated_value(T) )   // Omega(2,2)

      //!\return the value of the reduced dipole moment
      const CoeffType & delta_star() const;

    //! Friend base class so we can make implementation protected
    friend class SpeciesViscosityBase<KineticsTheoryViscosity<CoeffType,Interpolator>,CoeffType>;

    private:

      template <typename StateType>
      ANTIOCH_AUTO(StateType)
      op_impl(const StateType &T) const
      ANTIOCH_AUTOFUNC(StateType,  this->viscosity(T)  )

      void reset_coeffs_impl( const std::vector<CoeffType>& coeffs );

      void print_impl(std::ostream& os) const;

      void build_interpolation();

    //! Extrapolate to input maximum temperature, given in [K]
    /*!
     * The underlying collision integrals are interpolated for a given temperature range.
     * If the viscosity is to be evaluated outside that range, an error will occur.
     * This method will reconstruct the interpolation table, but use a linear extrapolation
     * from the max in the existing table to the input maximum temperature.
     */
    template <typename StateType>
    void extrapolate_max_temp_impl(const StateType& Tmax);

        //! building the spline
      void build_spline(const StockmayerPotential<CoeffType> & surface);

      /*! never ever use it*/
      KineticsTheoryViscosity();


      LennardJonesPotential<CoeffType> _LJ;
      CoeffType                        _dipole_moment;
      CoeffType                        _mass;

      Interpolator    _interp;
      CoeffType _delta_star;

        // a = 5 / 16 * sqrt( Boltzmann_constant * mass / pi) / ( sigma * sigma )
      CoeffType _a;

  };

  template <typename CoeffType, typename Interpolator>
  inline
  KineticsTheoryViscosity<CoeffType,Interpolator>::KineticsTheoryViscosity(const CoeffType & LJ_depth, const CoeffType & LJ_diameter,
                                                              const CoeffType & dipole_moment, const CoeffType & mass)
    : SpeciesViscosityBase<KineticsTheoryViscosity<CoeffType,Interpolator>,CoeffType>(),
    _LJ(LJ_depth,LJ_diameter),
    _dipole_moment(dipole_moment),
    _mass(mass),
    _delta_star(CoeffType(1e-7L) * ant_pow(Constants::light_celerity<CoeffType>(),2) * // * 1/(4*pi * eps_0) = 10^-7 * c^2
                ant_pow(_dipole_moment * Units<CoeffType>("D").get_SI_factor(),2) /
                ( _LJ.depth() * Constants::Boltzmann_constant<CoeffType>() * 2 * ant_pow(_LJ.diameter() * Units<CoeffType>("ang").get_SI_factor(),3) )),
    _a(CoeffType(0.3125e-14L) * ant_sqrt(CoeffType(1e28) * Constants::Boltzmann_constant<CoeffType>() * _mass / Constants::pi<CoeffType>())
                / (ant_pow(_LJ.diameter() * Units<CoeffType>("ang").get_SI_factor(),2))
       ) /* 5 / 16 * sqrt(pi * Boltzmann constant/pi) / (sigma^2)
            ~ 10^-14 float can't take 10^-28 in sqrt*/
  {
     this->build_interpolation();
     return;
  }

  template <typename CoeffType, typename Interpolator>
  inline
  KineticsTheoryViscosity<CoeffType,Interpolator>::KineticsTheoryViscosity(const std::vector<CoeffType> & coeffs):
#ifndef NDEBUG
    SpeciesViscosityBase<KineticsTheoryViscosity<CoeffType,Interpolator>,CoeffType>(),
        _LJ(-1,-1),
        _dipole_moment(-1),
        _mass(-1),
        _delta_star(-1),
       _a(0.L)
#else
    SpeciesViscosityBase<KineticsTheoryViscosity<CoeffType,Interpolator>,CoeffType>(),
        _LJ(coeffs[0],coeffs[1]),
        _dipole_moment(coeffs[2]),
        _mass(coeffs[3]),
        _delta_star(CoeffType(1e-7L) * ant_pow(Constants::light_celerity<CoeffType>(),2) * // * 1/(4*pi * eps_0) = 10^-7 * c^2
                    ant_pow(_dipole_moment * Units<CoeffType>("D").get_SI_factor(),2) /
                     ( _LJ.depth() * Constants::Boltzmann_constant<CoeffType>() * 2 * ant_pow(_LJ.diameter() * Units<CoeffType>("ang").get_SI_factor(),3) )),
        _a(CoeffType(0.3125e-14L) * ant_sqrt(CoeffType(1e28) * Constants::Boltzmann_constant<CoeffType>() * _mass / Constants::pi<CoeffType>())
                / (ant_pow(_LJ.diameter() * Units<CoeffType>("ang").get_SI_factor(),2))
          ) /* 5 / 16 * sqrt(pi * Boltzmann constant/pi) / (sigma^2)
                ~ 10^-14 float can't take 10^-28 in sqrt*/
#endif
  {
#ifndef NDEBUG
        antioch_assert_equal_to(coeffs.size(),4);

        this->reset_coeffs(coeffs[0],coeffs[1],coeffs[2],coeffs[3]);
#endif

     this->build_interpolation();
     return;
  }

  template <typename CoeffType, typename Interpolator>
  template <typename StateType>
  inline
  void KineticsTheoryViscosity<CoeffType,Interpolator>::extrapolate_max_temp_impl(const StateType & Tmax)
  {

     StockmayerPotential<CoeffType> surface;
    // Stockmayer is where the test is performed
     surface.extrapolate_to(Tmax/_LJ.depth());
     build_spline(surface);
  }

  template <typename CoeffType, typename Interpolator>
  inline
  void KineticsTheoryViscosity<CoeffType,Interpolator>::build_interpolation()
  {
     build_spline(StockmayerPotential<CoeffType>());
  }


  template <typename CoeffType, typename Interpolator>
  inline
  void KineticsTheoryViscosity<CoeffType,Interpolator>::build_spline(const StockmayerPotential<CoeffType> & surface)
  {
     std::vector<CoeffType> interp_surf(surface.log_temperature().size(),0);
     std::vector<CoeffType> rescaled_temp(surface.log_temperature().size(),0);
     for(unsigned int iT = 0; iT < surface.log_temperature().size(); iT++)
     {
        Interpolator spline(surface.delta(),surface.omega_2_2()[iT]);
        interp_surf[iT] = ant_sqrt(surface.temperature()[iT] * _LJ.depth()) /
                                spline.interpolated_value(_delta_star); // splining sqrt(T) / Omega<(2,2)>(log(T*))
        rescaled_temp[iT] = surface.temperature()[iT] * _LJ.depth();
     }

     _interp.spline_delete();
     _interp.spline_init(rescaled_temp,interp_surf); // T, sqrt(T)/Omega<(2,2)>(log(T*))
  }

  template <typename CoeffType, typename Interpolator>
  inline
  void KineticsTheoryViscosity<CoeffType,Interpolator>::reset_coeffs( const CoeffType & LJ_depth, const CoeffType & LJ_dia, const CoeffType & dipole_moment, const CoeffType & mass )
  {
//redefining parameters
     _LJ.reset_coeffs(LJ_depth,LJ_dia);
     _dipole_moment = dipole_moment;
     _mass = mass;
     _delta_star = CoeffType(1e-7L) * ant_pow(Constants::light_celerity<CoeffType>(),2) * // * 1/(4*pi * eps_0) = 10^-7 * c^2
                    ant_pow(_dipole_moment * Units<CoeffType>("D").get_SI_factor(),2) /
                     ( _LJ.depth() * Constants::Boltzmann_constant<CoeffType>() * 2 * ant_pow(_LJ.diameter() * Units<CoeffType>("ang").get_SI_factor(),3) );
           /* 5 / 16 * sqrt(pi * Boltzmann constant/pi) / (sigma^2)
                ~ 10^-14 float can't take 10^-28 in sqrt*/
     _a = CoeffType(0.3125e-14L) * ant_sqrt(CoeffType(1e28) * Constants::Boltzmann_constant<CoeffType>() * _mass / Constants::pi<CoeffType>())
                / (ant_pow(_LJ.diameter() * Units<CoeffType>("ang").get_SI_factor(),2));


//redefining collision integral
    this->build_interpolation();
  }

  template <typename CoeffType, typename Interpolator>
  inline
  void KineticsTheoryViscosity<CoeffType,Interpolator>::reset_coeffs_impl( const std::vector<CoeffType>& coeffs )
  {
     antioch_assert_equal_to(coeffs.size(),4);

     this->reset_coeffs(coeffs[0],coeffs[1],coeffs[2],coeffs[3]);
  }

  template <typename CoeffType, typename Interpolator>
  template <typename StateType>
  inline
  StateType KineticsTheoryViscosity<CoeffType,Interpolator>::compute_viscosity_and_derivative( const StateType& T, StateType & viscosity, StateType & dviscosity_dT ) const
  {
      // viscosity     = _a * spline(T)
      // dviscosity_dT = _a * dspline_dT(T)
     viscosity = this->viscosity(T);
     dviscosity_dT = _a * _interp.dinterp_dx(T);
     return;
  }

  template <typename CoeffType, typename Interpolator>
  inline
  void KineticsTheoryViscosity<CoeffType,Interpolator>::print_impl(std::ostream& os) const
  {
     os << "Pure species viscosity:\n"
        << "5/16 * sqrt(pi * " << _mass << " * kb * T) / ( pi * " << _LJ.depth() << "^2 * <Omega(2,2)*> )\n"
        << "T* = T / " << _LJ.depth() << "\n"
        << "delta* = 1/2 * " << _dipole_moment << "^2 / ( " << _LJ.depth() << " * kb * " << _LJ.diameter() << "^3 )\n"
        << "[factor a = " << _a << "]";
  }

  template <typename CoeffType, typename Interpolator>
  inline
  const CoeffType & KineticsTheoryViscosity<CoeffType,Interpolator>::delta_star() const
  {
     return _delta_star;
  }

} // end namespace Antioch

#endif //ANTIOCH_KINETICS_THEORY_VISCOSITY_H

#endif // ANTIOCH_HAVE_GSL
