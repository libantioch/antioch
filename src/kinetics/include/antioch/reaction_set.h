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

#ifndef ANTIOCH_REACTION_SET_H
#define ANTIOCH_REACTION_SET_H

// Antioch
#include "antioch/chemical_mixture.h"
#include "antioch/reaction.h"
#include "antioch/elementary_reaction.h"
#include "antioch/duplicate_reaction.h"
#include "antioch/threebody_reaction.h"
#include "antioch/falloff_reaction.h"
#include "antioch/lindemann_falloff.h"
#include "antioch/troe_falloff.h"

// C++
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

namespace Antioch
{
  /*!
   * This class encapsulates all the reaction mechanisms considered in a
   * chemical nonequilibrium simulation.
   */
  template<typename CoeffType=double>
  class ReactionSet
  {

  public:

    //! Constructor.
    ReactionSet( const ChemicalMixture<CoeffType>& chem_mixture );

    ~ReactionSet();

    //! \returns the number of species.
    unsigned int n_species() const;
     
    //! \returns the number of reactions.
    unsigned int n_reactions() const;
     
    //! Add a reaction to the system.
    void add_reaction(Reaction<CoeffType>* reaction);

    //! \returns a constant reference to reaction \p r.
    const Reaction<CoeffType>& reaction(const unsigned int r) const;

    const ChemicalMixture<CoeffType>& chemical_mixture() const;

    //! Compute the rates of progress for each reaction
    template <typename StateType, typename VectorStateType, typename VectorReactionsType>
    void compute_reaction_rates( const StateType& T,
                                 const VectorStateType& molar_densities,
                                 const VectorStateType& h_RT_minus_s_R,
                                 VectorReactionsType& net_reaction_rates ) const;

    //! Compute the rates of progress and derivatives for each reaction
    template <typename StateType, typename VectorStateType, typename VectorReactionsType, typename MatrixReactionsType>
    void compute_reaction_rates_and_derivs( const StateType& T,
                                            const VectorStateType& molar_densities,
                                            const VectorStateType& h_RT_minus_s_R,
                                            const VectorStateType& dh_RT_minus_s_R_dT,
                                            VectorReactionsType& net_reaction_rates,
                                            VectorReactionsType& dnet_rate_dT,
                                            MatrixReactionsType& dnet_rate_dX_s ) const;

    //!
    void update_particle_flux_chemistry();

    //!
    template <typename StateType, typename VectorStateType>
    void print_chemical_scheme( std::ostream& output,
                                const StateType& T,
                                const VectorStateType& molar_densities,
                                const VectorStateType& h_RT_minus_s_R ,
                                std::vector<VectorStateType> &lossMatrix,
                                std::vector<VectorStateType> &prodMatrix,
                                std::vector<VectorStateType> &netMatrix) const;

    //!
    template<typename StateType, typename VectorStateType>
    void get_reactive_scheme( const StateType& T,
                              const VectorStateType& molar_densities,
                              const VectorStateType& h_RT_minus_s_R,
                              VectorStateType& netRates,
                              VectorStateType& kfwdCoeff,
                              VectorStateType& kbkwdCoeff,
                              VectorStateType& kfwd,
                              VectorStateType& kbkwd,
                              VectorStateType& fwdC,
                              VectorStateType& bkwdC) const;


    //! Formatted print, by default to \p std::cout.
    void print( std::ostream& os = std::cout ) const;
     
    //! Formatted print.
    friend std::ostream& operator<<( std::ostream& os, const ReactionSet<CoeffType>& rset )
    {
      rset.print(os);
      return os;
    }

  private:
     
    ReactionSet();
     
    const ChemicalMixture<CoeffType>& _chem_mixture;

    std::vector<Reaction<CoeffType>* > _reactions;

    //! Scaling for equilibrium constant
    const CoeffType _P0_R;

  };
  
  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType>
  inline
  unsigned int ReactionSet<CoeffType>::n_species() const
  {
    return _chem_mixture.n_species();
  }

  template<typename CoeffType>
  inline
  unsigned int ReactionSet<CoeffType>::n_reactions() const
  {
    return _reactions.size();
  }

  template<typename CoeffType>
  inline
  void ReactionSet<CoeffType>::add_reaction(Reaction<CoeffType>* reaction)
  {
    _reactions.push_back(reaction);
    
    // and make sure it is initialized!
    _reactions.back()->initialize();

    return;
  }
  
  template<typename CoeffType>
  inline
  const Reaction<CoeffType>& ReactionSet<CoeffType>::reaction(const unsigned int r) const      
  {
    antioch_assert_less(r, this->n_reactions());
    return *_reactions[r];
  }

  template<typename CoeffType>
  inline
  const ChemicalMixture<CoeffType>& ReactionSet<CoeffType>::chemical_mixture() const
  {
    return _chem_mixture;
  }


  template<typename CoeffType>
  inline
  ReactionSet<CoeffType>::ReactionSet( const ChemicalMixture<CoeffType>& chem_mixture )
    : _chem_mixture(chem_mixture),
      _P0_R(1.0e5/Constants::R_universal<CoeffType>())
  {
    return;
  }


  template<typename CoeffType>
  inline
  ReactionSet<CoeffType>::~ReactionSet()
  {
    for(unsigned int ir = 1; ir < _reactions.size(); ir++)delete _reactions[ir];
    return;
  }
  

  template<typename CoeffType>
  inline
  void ReactionSet<CoeffType>::update_particle_flux_chemistry()
  {
      for(unsigned int ir = 0; ir < _reactions.size(); ir++)
      {
         _reactions[ir]->update_particles_flux();
      }
  }

  template<typename CoeffType>
  template<typename StateType, typename VectorStateType, typename VectorReactionsType>
  inline
  void ReactionSet<CoeffType>::compute_reaction_rates ( const StateType& T,
                                                        const VectorStateType& molar_densities,
                                                        const VectorStateType& h_RT_minus_s_R,
                                                        VectorReactionsType& net_reaction_rates ) const
  {
    antioch_assert_equal_to( net_reaction_rates.size(), this->n_reactions() );

    //!\todo Make these assertions vector-compatible
    // antioch_assert_greater(T, 0.0);

    antioch_assert_equal_to( molar_densities.size(), this->n_species() );
    antioch_assert_equal_to( h_RT_minus_s_R.size(), this->n_species() );

    // useful constants
    const StateType P0_RT = _P0_R/T; // used to transform equilibrium constant from pressure units

    // compute reaction forward rates & other reaction-sized arrays
    for (unsigned int rxn=0; rxn<this->n_reactions(); rxn++)
      {
        net_reaction_rates[rxn] = this->reaction(rxn).compute_rate_of_progress(molar_densities, T, P0_RT, h_RT_minus_s_R);
      }
    
    return;
  }

  template<typename CoeffType>
  template<typename StateType, typename VectorStateType, typename VectorReactionsType, typename MatrixReactionsType>
  inline
  void ReactionSet<CoeffType>::compute_reaction_rates_and_derivs( const StateType& T,
                                                                  const VectorStateType& molar_densities,
                                                                  const VectorStateType& h_RT_minus_s_R,
                                                                  const VectorStateType& dh_RT_minus_s_R_dT,
                                                                  VectorReactionsType& net_reaction_rates,
                                                                  VectorReactionsType& dnet_rate_dT,
                                                                  MatrixReactionsType& dnet_rate_dX_s ) const
  {
    antioch_assert_equal_to( net_reaction_rates.size(), this->n_reactions() );
    antioch_assert_equal_to( dnet_rate_dT.size(), this->n_reactions() );
    antioch_assert_equal_to( dnet_rate_dX_s.size(), this->n_reactions() );
#ifdef NDEBUG
#else
    for (unsigned int r=0; r < this->n_reactions(); r++)
      {
        antioch_assert_equal_to( dnet_rate_dX_s[r].size(), this->n_species() );
      }
#endif

    //!\todo Make these assertions vector-compatible
    // antioch_assert_greater(T, 0.0);

    antioch_assert_equal_to( molar_densities.size(), this->n_species() );
    antioch_assert_equal_to( h_RT_minus_s_R.size(), this->n_species() );

    // useful constants
    const StateType P0_RT = _P0_R/T; // used to transform equilibrium constant from pressure units

    // compute reaction forward rates & other reaction-sized arrays
    for (unsigned int rxn=0; rxn<this->n_reactions(); rxn++)
      {
        this->reaction(rxn).compute_rate_of_progress_and_derivatives( molar_densities, _chem_mixture, 
                                                                      T, P0_RT, h_RT_minus_s_R, dh_RT_minus_s_R_dT,
                                                                      net_reaction_rates[rxn], 
                                                                      dnet_rate_dT[rxn], 
                                                                      dnet_rate_dX_s[rxn] );
      }
    
    return;
  }

  template<typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void ReactionSet<CoeffType>::print_chemical_scheme( std::ostream& output,
                                                      const StateType& T,
                                                      const VectorStateType& molar_densities,
                                                      const VectorStateType& h_RT_minus_s_R,
                                                      std::vector<VectorStateType>& lossMatrix,
                                                      std::vector<VectorStateType>& prodMatrix,
                                                      std::vector<VectorStateType>& netMatrix ) const
  {

    //filling matrixes
    VectorStateType netRate,kfwdCoeff,kbkwdCoeff,kfwd,kbkwd,fwdC,bkwdC;

    //getting reaction infos
    get_reactive_scheme(T,molar_densities,h_RT_minus_s_R,netRate,kfwdCoeff,kbkwdCoeff,kfwd,kbkwd,fwdC,bkwdC);

    lossMatrix.resize(this->n_species());
    prodMatrix.resize(this->n_species());
    netMatrix.resize(this->n_species());

    for(unsigned int s = 0; s < this->n_species(); s++)
      {
        lossMatrix[s].resize(this->n_reactions(),0.);
        prodMatrix[s].resize(this->n_reactions(),0.);
        netMatrix[s].resize(this->n_reactions(),0.);
      }
    

    //filling matrixes
    for(unsigned int rxn = 0; rxn < this->n_reactions(); rxn++)
      {
        const Reaction<CoeffType>& reaction = this->reaction(rxn);
        for (unsigned int r=0; r<reaction.n_reactants(); r++)
          {
            lossMatrix[reaction.reactant_id(r)][rxn] += - static_cast<CoeffType>(reaction.reactant_stoichiometric_coefficient(r)) * kfwd[rxn];
            prodMatrix[reaction.reactant_id(r)][rxn] +=   static_cast<CoeffType>(reaction.reactant_stoichiometric_coefficient(r)) * kbkwd[rxn];
            netMatrix[reaction.reactant_id(r)][rxn]  +=   lossMatrix[reaction.reactant_id(r)][rxn] + prodMatrix[reaction.reactant_id(r)][rxn];
          }
        for (unsigned int p=0; p<reaction.n_products(); p++)
          {
            lossMatrix[reaction.product_id(p)][rxn] += - static_cast<CoeffType>(reaction.product_stoichiometric_coefficient(p)) * kbkwd[rxn];
            prodMatrix[reaction.product_id(p)][rxn] +=   static_cast<CoeffType>(reaction.product_stoichiometric_coefficient(p)) * kfwd[rxn];
            netMatrix[reaction.product_id(p)][rxn]  +=   lossMatrix[reaction.product_id(p)][rxn] + prodMatrix[reaction.product_id(p)][rxn];
          }
      }

    //explanation of header
    output << "# molar units considered for those budgets (kmol, m3, s)" << std::endl;
    output << "# formatted as follow:" << std::endl;
    output << "# rows: species, columns: reactions" << std::endl;
    output << "#" << std::endl;
    output << "#         |                        equation" << std::endl;
    output << "#         |  rate coefficient forward , rate coefficient backward" << std::endl;
    output << "#         |  forward concentrations   , backward concentrations" << std::endl;
    output << "#         |  rate forward             , rate backward" << std::endl;
    output << "# ------------------------------------------------------------------" << std::endl;
    output << "# species |  production term          ,  loss term" << std::endl;
    output << "# species |  net production/loss term" << std::endl;
    output << "#" << std::endl << std::endl;


    //header
    // // equation
    output << std::setw(10) << " ";
    for(unsigned int rxn = 0; rxn < this->n_reactions(); rxn++)
      {
        output << std::setw(31) << this->_reactions[rxn]->equation();
      }
    output << std::endl;
    // // coeffs
    output << std::setw(10) << " ";
    for(unsigned int rxn = 0; rxn < this->n_reactions(); rxn++)
      {
        output << std::setw(15) << kfwdCoeff[rxn] << "," << std::setw(15) << kbkwdCoeff[rxn];
      }
    output << std::endl;
    // // conc
    output << std::setw(10) << " ";
    for(unsigned int rxn = 0; rxn < this->n_reactions(); rxn++)
      {
        output << std::setw(15) << fwdC[rxn] << "," << std::setw(15) << bkwdC[rxn];
      }
    output << std::endl;
    // // rates
    output << std::setw(10) << " ";
    for(unsigned int rxn = 0; rxn < this->n_reactions(); rxn++)
      {
        output << std::setw(15) << kfwd[rxn] << "," << std::setw(15) << kbkwd[rxn];
      }
    output << std::endl;
    output << "----------------------------" << std::endl;

    //species

    for(unsigned int isp = 0; isp < this->n_species(); isp++)
      {
        output << std::setw(10) << _chem_mixture.chemical_species()[isp]->species();
        for(unsigned int rxn = 0; rxn < this->n_reactions(); rxn++)
          {
            output << std::setw(15) << std::scientific << std::setprecision(6) << prodMatrix[isp][rxn] << "," 
                   << std::setw(15) << std::scientific << std::setprecision(6) << lossMatrix[isp][rxn];
          }
        output << std::endl;
        // // net
        output << std::setw(10) << "";
        for(unsigned int rxn = 0; rxn < this->n_reactions(); rxn++)
          {
            output << std::setw(31) << std::scientific << std::setprecision(16) << netMatrix[isp][rxn];
          }
        output << std::endl;
      }
    output << std::endl;
    output << std::endl;
    output << std::endl;

    output << "# production table" << std::endl << std::endl;
    //header
    // // equation
    output << std::setw(10) << " ";
    for(unsigned int rxn = 0; rxn < this->n_reactions(); rxn++)
      {
        output << std::setw(31) << this->_reactions[rxn]->equation();
      }
    output << std::setw(31) << "Total" << std::endl;
    for(unsigned int isp = 0; isp < this->n_species(); isp++)
      {
        output << std::setw(10) << _chem_mixture.chemical_species()[isp]->species();
        CoeffType rowTotal(0.);
        for(unsigned int rxn = 0; rxn < this->n_reactions(); rxn++)
          {
            output << std::setw(31) << std::scientific << std::setprecision(16) << prodMatrix[isp][rxn];
            rowTotal += prodMatrix[isp][rxn];
          }
        output << std::setw(31) << std::scientific << std::setprecision(16) << rowTotal << std::endl;
      }
    output << std::endl << std::endl;

    output << "# loss table" << std::endl << std::endl;
    //header
    // // equation
    output << std::setw(10) << " ";
    for(unsigned int rxn = 0; rxn < this->n_reactions(); rxn++)
      {
        output << std::setw(31) << this->_reactions[rxn]->equation();
      }
    output << std::setw(31) << "Total" << std::endl;
    for(unsigned int isp = 0; isp < this->n_species(); isp++)
      {
        output << std::setw(10) << _chem_mixture.chemical_species()[isp]->species();
        CoeffType rowTotal(0.);
        for(unsigned int rxn = 0; rxn < this->n_reactions(); rxn++)
          {
            output << std::setw(31) << std::scientific << std::setprecision(16) << lossMatrix[isp][rxn];
            rowTotal += lossMatrix[isp][rxn];
          }
        output << std::setw(31) << std::scientific << std::setprecision(16) << rowTotal << std::endl;
      }
    output << std::endl << std::endl;

    output << "# net table" << std::endl << std::endl;
    //header
    // // equation
    output << std::setw(10) << " ";
    for(unsigned int rxn = 0; rxn < this->n_reactions(); rxn++)
      {
        output << std::setw(31) << this->_reactions[rxn]->equation();
      }
    output << std::setw(31) << "Total" << std::endl;

    CoeffType columnTotal(0.);
    std::vector<CoeffType> columnSum;
    columnSum.resize(this->n_reactions(),0.);
    CoeffType netTotalRow(0.);

    for(unsigned int isp = 0; isp < this->n_species(); isp++)
      {
        output << std::setw(10) << _chem_mixture.chemical_species()[isp]->species();
        CoeffType rowTotal(0.);
        for(unsigned int rxn = 0; rxn < this->n_reactions(); rxn++)
          {
            output << std::setw(31) << std::scientific << std::setprecision(16) << netMatrix[isp][rxn];
            rowTotal += netMatrix[isp][rxn];
            columnSum[rxn] += netMatrix[isp][rxn];
          }
        output << std::setw(31) << std::scientific << std::setprecision(16) << rowTotal << std::endl;
        netTotalRow += rowTotal;
      }
    output << std::setw(10) << "Total";
    for(unsigned int rxn = 0; rxn < this->n_reactions(); rxn++)
      {
        output << std::setw(31) << std::scientific << std::setprecision(16) << columnSum[rxn];
        columnTotal += columnSum[rxn];
      }
    output << std::endl << std::endl;
    output << "sum of row sums:    " << std::scientific << std::setprecision(16) << netTotalRow << std::endl;
    output << "sum of column sums: " << std::scientific << std::setprecision(16) << columnTotal << std::endl << std::endl;
    output << "# net table, mass here (kg, m3, s)" << std::endl << std::endl;
    //header
    // // equation
    output << std::setw(10) << " ";
    for(unsigned int rxn = 0; rxn < this->n_reactions(); rxn++)
      {
        output << std::setw(31) << this->_reactions[rxn]->equation();
      }
    output << std::setw(31) << "Total" << std::endl;

    columnTotal = 0.;
    columnSum.resize(this->n_reactions(),0.);
    netTotalRow = 0.;

    for(unsigned int isp = 0; isp < this->n_species(); isp++)
      {
        output << std::setw(10) << _chem_mixture.chemical_species()[isp]->species();
        CoeffType rowTotal(0.);
        for(unsigned int rxn = 0; rxn < this->n_reactions(); rxn++)
          {
            output << std::setw(31) << std::scientific << std::setprecision(16) << netMatrix[isp][rxn];
            rowTotal += netMatrix[isp][rxn] * _chem_mixture.M(isp);
            columnSum[rxn] += netMatrix[isp][rxn] * _chem_mixture.M(isp);
          }
        output << std::setw(31) << std::scientific << std::setprecision(16) << rowTotal << std::endl;
        netTotalRow += rowTotal;
      }
    output << std::setw(10) << "Total";
    for(unsigned int rxn = 0; rxn < this->n_reactions(); rxn++)
      {
        output << std::setw(31) << std::scientific << std::setprecision(16) << columnSum[rxn];
        columnTotal += columnSum[rxn];
      }
    output << std::endl << std::endl;

    output << std::endl << std::endl;
    output << "# Mixture informations" << std::endl;
    output << "# species / concentration / molar mass" << std::endl;
    for(unsigned int isp = 0; isp < this->n_species(); isp++)
      {
        output << std::setw(10) << _chem_mixture.chemical_species()[isp]->species()
               << std::setw(31) << molar_densities[isp]
               << std::setw(31) << _chem_mixture.M(isp) << std::endl;
      }

    return;
  }

  template<typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void ReactionSet<CoeffType>::get_reactive_scheme( const StateType& T,
                                                    const VectorStateType& molar_densities,
                                                    const VectorStateType& h_RT_minus_s_R,
                                                    VectorStateType& netRates,
                                                    VectorStateType& kfwdCoeff,
                                                    VectorStateType& kbkwdCoeff,
                                                    VectorStateType& kfwd,
                                                    VectorStateType& kbkwd,
                                                    VectorStateType& fwdC,
                                                    VectorStateType& bkwdC) const
  {
    //!\todo Make these assertions vector-compatible
    // antioch_assert_greater(T, 0.0);

    antioch_assert_equal_to( molar_densities.size(), this->n_species() );
    antioch_assert_equal_to( h_RT_minus_s_R.size(), this->n_species() );

    // useful constants
    const StateType P0_RT = _P0_R/T; // used to transform equilibrium constant from pressure units

    netRates.resize(this->n_reactions(),0.);
    kfwdCoeff.resize(this->n_reactions(),0.);
    kfwd.resize(this->n_reactions(),0.);
    kbkwdCoeff.resize(this->n_reactions(),0.);
    kbkwd.resize(this->n_reactions(),0.);
    fwdC.resize(this->n_reactions(),1.);
    bkwdC.resize(this->n_reactions(),1.);

    // compute reaction forward rates & other reaction-sized arrays
    for (unsigned int rxn=0; rxn<this->n_reactions(); rxn++)
      {
        const Reaction<CoeffType>& reaction = this->reaction(rxn);
        kfwdCoeff[rxn] = reaction.compute_forward_rate_coefficient(molar_densities,T);
        kfwd[rxn] = kfwdCoeff[rxn];
        StateType keq = reaction.equilibrium_constant( P0_RT, h_RT_minus_s_R );

        kbkwdCoeff[rxn] = kfwdCoeff[rxn]/keq;
        kbkwd[rxn] = kbkwdCoeff[rxn];

        for (unsigned int r=0; r<reaction.n_reactants(); r++)
          {
            fwdC[rxn] *= pow( molar_densities[reaction.reactant_id(r)],
                              static_cast<int>(reaction.reactant_stoichiometric_coefficient(r)) );
          }
        kfwd[rxn] *= fwdC[rxn];
        for (unsigned int p=0; p<reaction.n_products(); p++)
          {
            bkwdC[rxn] *= pow( molar_densities[reaction.product_id(p)],
                               static_cast<int>(reaction.product_stoichiometric_coefficient(p)) );
          }
        kbkwd[rxn] *= bkwdC[rxn];

        netRates[rxn] = kfwd[rxn] - kbkwd[rxn];

      }

    return;
  }


  template<typename CoeffType>
  inline
  void ReactionSet<CoeffType>::print(std::ostream& os) const
  {
    os << "# Number of reactions: " << this->n_reactions() << "\n";

    for (unsigned int r=0; r < this->n_reactions(); r++)
      {
        os << "# " << r << '\n'
           << (*this->reaction(r)) << "\n";
      }

    return;
  }
 
} // end namespace Antioch

#endif // ANTIOCH_REACTION_SET_H
