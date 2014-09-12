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

#ifndef ANTIOCH_MOLECULAR_BINARY_DIFFUSION_H
#define ANTIOCH_MOLECULAR_BINARY_DIFFUSION_H

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/math_constants.h"
#include "antioch/cmath_shims.h"
#include "antioch/Stockmayer_potential.h"
#include "antioch/transport_species.h"

//C++

namespace Antioch{

  template <typename CoeffType>
  class MolecularBinaryDiffusion
  {
        public:
          MolecularBinaryDiffusion();
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

          //! \return molecular binary diffusion coefficient
          template <typename StateType, typename VectorStateType, typename SpeciesType>
          const 
          ANTIOCH_AUTO(StateType)
          operator()(const StateType & T, const VectorStateType & molar_concentrations, 
                     const VectorStateType & collision_interp,
                     const TransportSpecies<SpeciesType> & si, const TransportSpecies<SpeciesType> & sj) const
          ANTIOCH_AUTOFUNC(StateType, this->binary_diffusion(T,molar_concentrations,collision_interp,si,sj))

          //! \return molecular binary diffusion coefficient
          template <typename StateType, typename VectorStateType, typename SpeciesType>
          const StateType binary_diffusion(const StateType & T, const VectorStateType & molar_concentrations, 
                                           const VectorStateType & collision_interp,
                                           const TransportSpecies<SpeciesType> & si, const TransportSpecies<SpeciesType> & sj) const;


          DiffusionModel model() const {return _model;}


        private:

          template <typename SpeciesType>
          const CoeffType composed_xi(const TransportSpecies<SpeciesType> & si, const TransportSpecies<SpeciesType> & sj) const;

          const CoeffType _coefficient;

          const DiffusionModel _model;

  };

  template <typename CoeffType>
  inline
  MolecularBinaryDiffusion<CoeffType>::~MolecularBinaryDiffusion()
  {
     return;
  }

  template <typename CoeffType>
  inline
  MolecularBinaryDiffusion<CoeffType>::MolecularBinaryDiffusion():
/* cheating on the powers*/
      _coefficient(CoeffType(3e-3L * 4.L /16.L) * ant_sqrt(CoeffType(2.L) * (Constants::Boltzmann_constant<CoeffType>() * CoeffType(1e23L) ) /
                                                            ((Constants::Avogadro<CoeffType>() * CoeffType(1e-23)) * Constants::pi<CoeffType>())
                                                          )
                  ), // 4 is from reduced diameter => 1 / (1/2)^2, 10^20 for SI, ang -> m
      _model(BIMOLECULAR)
  {
     return;
  }

  template <typename CoeffType>
  template <typename SpeciesType>
  inline
  const CoeffType MolecularBinaryDiffusion<CoeffType>::composed_xi(const TransportSpecies<SpeciesType> & si, const TransportSpecies<SpeciesType> & sj) const
  {
        const TransportSpecies<SpeciesType> &n = (si.polar())?sj:si;
        const TransportSpecies<SpeciesType> &p = (si.polar())?si:sj;

        CoeffType pol    =  n.polarizability() / ant_pow(n.LJ_diameter(),3); //ang^3 / ang^3
// mu/sqrt(4 pi eps_0 eps sig^3) = mu * sqrt(10^7 * c^2 / ( eps sig^3))
        CoeffType dipole = p.dipole_moment() * CoeffType(3.335641e-30L) * // Debye
                           ant_sqrt(Constants::light_celerity<CoeffType>() * Constants::light_celerity<CoeffType>() * CoeffType(1e23L) // 1/(4 * pi * eps_0) =  10^7 * c^2, to SI ang^3 -> m^3
                                     / ( p.LJ_depth() * Constants::Boltzmann_constant<CoeffType>() * ant_pow(p.LJ_diameter(),3) ) );

        return (CoeffType)(1.L) + (CoeffType)(0.25L) * pol * dipole * ant_sqrt(p.LJ_depth()/n.LJ_depth());
  }

  template <typename CoeffType>
  template <typename StateType, typename VectorStateType, typename SpeciesType>
  inline
  const StateType MolecularBinaryDiffusion<CoeffType>::binary_diffusion(const StateType & T, const VectorStateType & molar_concentrations,
                                                                        const VectorStateType & collision_interp,
                                                                        const TransportSpecies<SpeciesType> & si, const TransportSpecies<SpeciesType> & sj) const
  { 
     StateType xi  = (si.polar() == sj.polar())?1:this->composed_xi(si,sj);
     StateType red_mass = (si.species() == sj.species())?si.M():(si.M() * sj.M()) / (si.M() + sj.M()); // kg/mol

     StateType nTot(0);
     for(unsigned int s = 0; s < molar_concentrations.size(); s++)
     {
         nTot += molar_concentrations[s];
     }

// collision integral
     StateType omega_11(0);
     StateType T_star(T / (ant_sqrt(si.LJ_depth() * sj.LJ_depth()) * xi * xi));
     StateType Tmp(1);
     StateType log_T_star = ant_log(T_star);
// on log of reduced temperature
     for(unsigned int d = 0; d < collision_interp.size(); d++)
     {
        omega_11 += collision_interp[d] * Tmp;
        Tmp *= log_T_star;
     }

     return _coefficient * ant_sqrt(T / red_mass) / (nTot * ant_pow(si.LJ_diameter() + sj.LJ_diameter(),2) * ant_pow(xi,-1.L/3.L) * omega_11); 
  }

}

#endif
