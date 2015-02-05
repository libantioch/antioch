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

#include "antioch_config.h"
#ifdef ANTIOCH_HAVE_GSL // if we do not have it, we don't even define the stuff

#ifndef ANTIOCH_MOLECULAR_BINARY_DIFFUSION_H
#define ANTIOCH_MOLECULAR_BINARY_DIFFUSION_H


//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/math_constants.h"
#include "antioch/cmath_shims.h"
#include "antioch/Stockmayer_potential.h"
#include "antioch/transport_species.h"
#include "antioch/gsl_spliner.h"

//C++
#include <vector>



namespace Antioch{

  template <typename CoeffType, typename Interpolator = GSLSpliner>
  class MolecularBinaryDiffusion
  {
        typedef unsigned int Species;

        public:
          MolecularBinaryDiffusion(const TransportSpecies<CoeffType> & si, const TransportSpecies<CoeffType> & sj);
          MolecularBinaryDiffusion(const std::vector<TransportSpecies<CoeffType> >& species);
          ~MolecularBinaryDiffusion();

          /*!
           *
           * \TODO derivatives
           *
           * \f[
           *   D_{ij} = \frac{3}{16} \sqrt{\frac{2 \mathrm{k_B}^3 \mathrm{T}^3}{\pi m_{ij}}}\frac{1}{\mathrm{P} \sigma_{ij}^2\Omega^{(1,1)*}}
           *          = \frac{3\sqrt{2 \mathrm{k_B}^2}}{16\sqrt{\pi}\mathcal{N}_\text{Avo}} \frac{T^2}{\sum_s c_s m_{ij} \sigma_{ij}^2\Omega^{(1,1)*}}
           * \f]
           * with \f$T\f$ and \f$P\f$ the temperature and pressure respectively, \f$m_{ij}\f$ the reduced mass, given by
           * \f]
           *   \begin{array}{l@{,\quad}r}
           *      m_{ij} = \frac{m_jm_i}{m_i + m_j} & \text{ if } j \neq i \\
           *      m_{ij} = m_i                      & \text{ if } j = i    \\
           *    \end{array}
           * \f]
           * \f$\sigma_{ij}\f$ the reduced collision diameter and \f$\Omega^{(1,1)*}\f$ the collision integral.
           *  
           * Using the ideal gas approximation \f$\mathrm{P} = n \mathrm{RT}\f$ with \f$n\f$ the molar concentration
           * (SI units are \f$\mathtt{mol\,m^{-3}}\f$). Thus we obtain
           * \f[
           *   D_{ij} = \frac{3}{16} \frac{\sqrt{2\pi \mathrm{k_B}^3 \mathrm{T}^2\frac{1}{m_{ij}}}}{n_{\mathrm{tot}}\mathrm{R} \pi \sigma_{jk}^2\Omega^{(1,1)*}}
           * \f]
           *  with $\f$n_{\mathrm{tot}}\f$ the total concentration of the mixture:
           * \f[
           *   n_{\mathrm{tot}} = \sum_{s} n_s
           * \f]
           */

          void reset_coeffs(const std::vector<TransportSpecies<CoeffType> > & spec);
          void reset_coeffs(const TransportSpecies<CoeffType> & si, const TransportSpecies<CoeffType> & sj);

          //! \return molecular binary diffusion coefficient
          template <typename StateType>
          const 
          ANTIOCH_AUTO(StateType)
          operator()(const StateType & T, const StateType & cTot) const
          ANTIOCH_AUTOFUNC(StateType, this->binary_diffusion(T, cTot))

          template <typename StateType>
          const
          ANTIOCH_AUTO(StateType)
             Stockmayer(const StateType & T)
          ANTIOCH_AUTOFUNC(StateType, ant_sqrt(T) / _interp.interpolated_value(T) )

          //! \return molecular binary diffusion coefficient
          template <typename StateType>
          const 
          ANTIOCH_AUTO(StateType) 
                binary_diffusion(const StateType & T, const StateType & cTot) const
          ANTIOCH_AUTOFUNC(StateType,_coefficient * _interp.interpolated_value(T) / cTot )

           //!
          const CoeffType & reduced_dipole_moment() const;

           //!
          const CoeffType & reduced_mass() const;

           //!
          const CoeffType & reduced_LJ_diameter() const;

           //!
          const CoeffType & reduced_LJ_depth() const;

           //!
          const CoeffType & xi() const;

          void print(std::ostream & out = std::cout) const;

          friend std::ostream & operator<< (std::ostream & out, const MolecularBinaryDiffusion<CoeffType,Interpolator> & bimol_diff)
          {
            bimol_diff.print(out);
            return out;
          }
        

        private:

          MolecularBinaryDiffusion();

          // convenient method to avoid code duplication
          void set_coeffs(const std::vector<TransportSpecies<CoeffType> > & spec);

          template <typename SpeciesType>
          const CoeffType composed_xi(const TransportSpecies<SpeciesType> & si, const TransportSpecies<SpeciesType> & sj);

          void build_interpolation();


          Interpolator _interp;

// parameter for the couple
          Species   _i,_j;
          CoeffType _reduced_mass;
          CoeffType _xi; // simplify calculations
          CoeffType _reduced_LJ_diameter;
          CoeffType _reduced_LJ_depth;
          CoeffType _reduced_dipole_moment;

          CoeffType _coefficient;

  };

  template <typename CoeffType, typename Interpolator>
  inline
  MolecularBinaryDiffusion<CoeffType,Interpolator>::~MolecularBinaryDiffusion()
  {
     return;
  }

  template <typename CoeffType, typename Interpolator>
  inline
  MolecularBinaryDiffusion<CoeffType,Interpolator>::MolecularBinaryDiffusion(const TransportSpecies<CoeffType> & si, const TransportSpecies<CoeffType> & sj):
     _i(si.species()),
     _j(sj.species()),
     _reduced_mass((si.species() == sj.species())?CoeffType(0.5L) * si.M():(si.M() * sj.M()) / (si.M() + sj.M())), // kg/mol
     _xi((si.polar() == sj.polar())?1L:this->composed_xi(si,sj)),
     _reduced_LJ_diameter(CoeffType(0.5L) * (si.LJ_diameter() + sj.LJ_diameter()) * Units<CoeffType>("ang").get_SI_factor() * ant_pow(_xi,-1.L/6.L)), // 1/2 * (sigma_1 + sigma_2) * xi^(-1/6)
     _reduced_LJ_depth(ant_sqrt(si.LJ_depth() * sj.LJ_depth())  * _xi * _xi), // sqrt(eps_1 * eps_2) * xi^2
     _reduced_dipole_moment((CoeffType(1e-7L) * ant_pow(Constants::light_celerity<CoeffType>(),2)) * 
                              (si.dipole_moment() * sj.dipole_moment() * ant_pow(Units<CoeffType>("D").get_SI_factor(),2 )) / 
                            (CoeffType(2.L) * _reduced_LJ_depth * Constants::Boltzmann_constant<CoeffType>() * ant_pow(_reduced_LJ_diameter,3) )), // mu^2 / (2*eps*sigma^3)
/*  = ~ 7.16 10^-25, 
        float can't take it, 
        we cheat on the kb / nAvo division that makes float cry

        3/16 * sqrt ( 2 * kb / (Navo * pi) )
*/
      _coefficient(CoeffType(0.1875e-25L) * ant_sqrt(CoeffType(2.L) * Constants::Boltzmann_constant<CoeffType>() * CoeffType(1e25L)  /
                                                            (Constants::Avogadro<CoeffType>() * CoeffType(1e-25L) * Constants::pi<CoeffType>())
                                                          )
                                         / ( sqrt(_reduced_mass) * (_reduced_LJ_diameter * _reduced_LJ_diameter) )
                  )
  {
     
     this->build_interpolation();

     return;
  }

  template <typename CoeffType, typename Interpolator>
  inline
  MolecularBinaryDiffusion<CoeffType,Interpolator>::MolecularBinaryDiffusion(const std::vector<TransportSpecies<CoeffType> >& species):
#ifndef NDEBUG
     _i(0),
     _j(0),
     _reduced_mass(-1),
     _reduced_dipole_moment(-1),
     _xi(-1),
     _reduced_LJ_diameter(-1),
     _reduced_LJ_depth(-1),
      _coefficient(0.)
#else
     _i(species[0].species()),
     _j(species[1].species()),
     _reduced_mass((species[0].species() == species[1].species())?CoeffType(0.5L) * species[0].M():(species[0].M() * species[1].M()) / (species[0].M() + species[1].M())), // kg/mol
     _xi((species[0].polar() == species[1].polar())?1L:this->composed_xi(species[0],species[1])),
     _reduced_LJ_diameter(CoeffType(0.5L) * (species[0].LJ_diameter() + species[1].LJ_diameter()) * Units<CoeffType>("ang").get_SI_factor() * ant_pow(_xi,-1.L/6.L)), // 1/2 * (sigma_1 + sigma_2) * xi^(-1/6)
     _reduced_LJ_depth(ant_sqrt(species[0].LJ_depth() * species[1].LJ_depth())  *_xi * _xi), // sqrt(eps_1 * eps_2) * xi^2
     _reduced_dipole_moment((CoeffType(1e-7L) * ant_pow(Constants::light_celerity<CoeffType>(),2)) * 
                              (species[0].dipole_moment() * species[1].dipole_moment() * ant_pow(Units<CoeffType>("D").get_SI_factor(),2 )) / 
                            (CoeffType(2.L) * _reduced_LJ_depth * Constants::Boltzmann_constant<CoeffType>() * ant_pow(_reduced_LJ_diameter,3) )), // mu^2 / (2*eps*sigma^3)
      _coefficient(CoeffType(0.1875e-25L) * ant_sqrt(CoeffType(2.L) * Constants::Boltzmann_constant<CoeffType>() * CoeffType(1e25L)  /
                                                            (Constants::Avogadro<CoeffType>() * CoeffType(1e-25L) * Constants::pi<CoeffType>())
                                                          )
                                         / ( sqrt(_reduced_mass) * (_reduced_LJ_diameter * _reduced_LJ_diameter) )
                  )
#endif
  {
#ifndef NDEBUG
     antioch_assert_equal_to(species.size(),2);

     this->set_coeffs(species);
#endif

     this->build_interpolation();
  }

  template <typename CoeffType, typename Interpolator>
  inline
  void MolecularBinaryDiffusion<CoeffType,Interpolator>::reset_coeffs(const std::vector<TransportSpecies<CoeffType> >& species)
  {
     antioch_assert_equal_to(species.size(),2);

     this->reset_coeffs(species[0],species[1]);
  }

  template <typename CoeffType, typename Interpolator>
  inline
  void MolecularBinaryDiffusion<CoeffType,Interpolator>::reset_coeffs(const TransportSpecies<CoeffType> & si, const TransportSpecies<CoeffType> & sj)
  {
     _i                     = si.species();
     _j                     = sj.species();
     _reduced_mass          = (si.species() == sj.species())?si.M():(si.M() * sj.M()) / (si.M() + sj.M());
     _reduced_dipole_moment = ant_sqrt(si.dipole_moment() * sj.dipole_moment());
     _xi                    = (si.polar() == sj.polar())?1L:this->composed_xi(si,sj);
     _reduced_LJ_diameter   = CoeffType(0.5L) * (si.LJ_diameter() + sj.LJ_diameter()) * Units<CoeffType>("ang").get_SI_factor() * ant_pow(_xi,-1.L/6.L);
     _reduced_LJ_depth      = ant_sqrt(si.LJ_depth() * sj.LJ_depth()) *_xi * _xi;
     _reduced_dipole_moment = CoeffType(1e-7L) * ant_pow(Constants::light_celerity<CoeffType>(),2) * // * 1/(4*pi * eps_0) = 10^-7 * c^2 
                             (si.dipole_moment() * sj.dipole_moment()) * ant_pow(Units<CoeffType>("D").get_SI_factor(),2)
                                / ((CoeffType(2.L)*_reduced_LJ_depth * Constants::Boltzmann_constant<CoeffType>() * ant_pow(_reduced_LJ_diameter,3)));
/*  = ~ 7.16 10^-25, 
        float can't take it, 
        we cheat on the kb / nAvo division that makes float cry

        3/16 * sqrt ( 2 * kb / (Navo * pi) )
*/
      _coefficient = CoeffType(0.1875e-25L) * ant_sqrt(CoeffType(2.L) * Constants::Boltzmann_constant<CoeffType>() * CoeffType(1e25L)  /
                                                            (Constants::Avogadro<CoeffType>() * CoeffType(1e-25L) * Constants::pi<CoeffType>())
                                                          )
                                         / ( sqrt(_reduced_mass) * (_reduced_LJ_diameter * _reduced_LJ_diameter) );
  }

  template <typename CoeffType, typename Interpolator>
  inline
  void MolecularBinaryDiffusion<CoeffType,Interpolator>::set_coeffs(const std::vector<TransportSpecies<CoeffType> >& species)
  {
     antioch_assert_equal_to(species.size(),2);

     this->reset_coeffs(species[0],species[1]);
  }

  template <typename CoeffType, typename Interpolator>
  template <typename SpeciesType>
  inline
  const CoeffType MolecularBinaryDiffusion<CoeffType,Interpolator>::composed_xi(const TransportSpecies<SpeciesType> & si, const TransportSpecies<SpeciesType> & sj)
  {
        const TransportSpecies<SpeciesType> &n = (si.polar())?sj:si;
        const TransportSpecies<SpeciesType> &p = (si.polar())?si:sj;

        CoeffType pol    =  n.polarizability() / ant_pow(n.LJ_diameter(),3); //ang^3 / ang^3 -> cancel out
        CoeffType dipole = p.dipole_moment() * Units<CoeffType>("D").get_SI_factor() 
                           / ant_sqrt((CoeffType)4.L * Constants::pi<CoeffType>() * Constants::vacuum_permittivity<CoeffType>()
                                        *  p.LJ_depth() * ant_pow(p.LJ_diameter(),3) );

        return (CoeffType)(1.L) + (CoeffType)(0.25L) * pol * dipole * ant_sqrt(p.LJ_depth()/n.LJ_depth());
  }

  template <typename CoeffType, typename Interpolator>
  inline
  void MolecularBinaryDiffusion<CoeffType,Interpolator>::build_interpolation()
  {
     _interp.spline_delete();
     StockmayerPotential<CoeffType> surface;
     std::vector<CoeffType> interp_surf(surface.log_temperature().size(),0);
     std::vector<CoeffType> rescaled_temp(surface.log_temperature().size(),0);
     for(unsigned int iT = 0; iT < surface.temperature().size(); iT++)
     {
        Interpolator spline(surface.delta(),surface.omega_1_1()[iT]);
        interp_surf[iT] = ant_sqrt(surface.temperature()[iT] * _reduced_LJ_depth) / 
                                spline.interpolated_value(_reduced_dipole_moment); // splining sqrt(T) / Omega<(1,1)>(log(T*))
        rescaled_temp[iT] = surface.temperature()[iT] * _reduced_LJ_depth;
     }

     _interp.spline_init(rescaled_temp,interp_surf);
  }

  template <typename CoeffType, typename Interpolator>
  inline
  void MolecularBinaryDiffusion<CoeffType,Interpolator>::print(std::ostream & out) const
  {
     out << "Molecular binary diffusion of couple " << _i << "," << _j << ":\n"
         << "  reduced mass = "          << _reduced_mass          << "\n"
         << "  reduced dipole moment = " << _reduced_dipole_moment << "\n"
         << "  xi = "                    << _xi                    << "\n"
         << "  reduced LJ diameter = "   << _reduced_LJ_diameter   << "\n"
         << "  reduced LJ depth = "      << _reduced_LJ_depth 
         << std::endl;
  }

  template <typename CoeffType, typename Interpolator>
  inline
  const CoeffType & MolecularBinaryDiffusion<CoeffType,Interpolator>::reduced_dipole_moment() const
  {
     return _reduced_dipole_moment;
  }

  template <typename CoeffType, typename Interpolator>
  inline
  const CoeffType & MolecularBinaryDiffusion<CoeffType,Interpolator>::reduced_mass() const
  {
     return _reduced_mass;
  }

  template <typename CoeffType, typename Interpolator>
  inline
  const CoeffType & MolecularBinaryDiffusion<CoeffType,Interpolator>::reduced_LJ_diameter() const
  {
     return _reduced_LJ_diameter;
  }

  template <typename CoeffType, typename Interpolator>
  inline
  const CoeffType & MolecularBinaryDiffusion<CoeffType,Interpolator>::reduced_LJ_depth() const
  {
     return _reduced_LJ_depth;
  }

  template <typename CoeffType, typename Interpolator>
  inline
  const CoeffType & MolecularBinaryDiffusion<CoeffType,Interpolator>::xi() const
  {
     return _xi;
  }


}

#include "antioch/molecular_binary_diffusion_utils_decl.h"
#endif // ANTIOCH_BIMOL_DIFF

#endif // ifdef ANTIOCH_HAVE_GSL
