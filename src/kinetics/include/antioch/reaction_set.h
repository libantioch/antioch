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

#ifndef ANTIOCH_REACTION_SET_H
#define ANTIOCH_REACTION_SET_H

// Antioch
#include "antioch/chemical_mixture.h"
#include "antioch/reaction.h"
#include "antioch/kinetics_conditions.h"
#include "antioch/elementary_reaction.h"
#include "antioch/duplicate_reaction.h"
#include "antioch/threebody_reaction.h"
#include "antioch/falloff_reaction.h"
#include "antioch/falloff_threebody_reaction.h"
#include "antioch/lindemann_falloff.h"
#include "antioch/troe_falloff.h"
#include "antioch/string_utils.h"

// C++
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <limits>

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
    //
    // The ownership is transfered to the ReactionSet
    // object, they are deleted in the destructor.
    void add_reaction(Reaction<CoeffType>* reaction);

    //! remove a reaction from the system.
    //
    // The corresponding pointer is deleted
    void remove_reaction(unsigned int nr);

    //! \returns a constant reference to reaction \p r.
    const Reaction<CoeffType>& reaction(const unsigned int r) const;

    //! \returns a writeable reference to reaction \p r.
    Reaction<CoeffType>& reaction(const unsigned int r);

    //! \returns the index of a reaction given its id
    unsigned int reaction_by_id(const std::string & reaction_id) const;

    //! change a parameter of a reaction
    //
    // in charge of the human-to-antioch translation
    template <typename ParamType>
    void set_parameter_of_reaction(const std::string & reaction_id, const std::vector<std::string> & keywords, ParamType value);

    //! \return a parameter of a reaction
    //
    // in charge of the human-to-antioch translation
    CoeffType get_parameter_of_reaction(const std::string & reaction_id, const std::vector<std::string> & keywords) const;

    const ChemicalMixture<CoeffType>& chemical_mixture() const;

    //! Compute the rates of progress for each reaction
    template <typename StateType, typename VectorStateType, typename VectorReactionsType>
    void compute_reaction_rates( const KineticsConditions<StateType,VectorStateType>& conditions,
                                 const VectorStateType& molar_densities,
                                 const VectorStateType& h_RT_minus_s_R,
                                 VectorReactionsType& net_reaction_rates ) const;

    //! Compute the rates of progress and derivatives for each reaction
    template <typename StateType, typename VectorStateType, typename VectorReactionsType, typename MatrixReactionsType>
    void compute_reaction_rates_and_derivs( const KineticsConditions<StateType,VectorStateType>& conditions,
                                            const VectorStateType& molar_densities,
                                            const VectorStateType& h_RT_minus_s_R,
                                            const VectorStateType& dh_RT_minus_s_R_dT,
                                            VectorReactionsType& net_reaction_rates,
                                            VectorReactionsType& dnet_rate_dT,
                                            MatrixReactionsType& dnet_rate_dX_s ) const;

    //!
    template <typename StateType, typename VectorStateType>
    void print_chemical_scheme( std::ostream& output,
                                const KineticsConditions<StateType,VectorStateType>& conditions,
                                const VectorStateType& molar_densities,
                                const VectorStateType& h_RT_minus_s_R ,
                                std::vector<VectorStateType> &lossMatrix,
                                std::vector<VectorStateType> &prodMatrix,
                                std::vector<VectorStateType> &netMatrix) const;

    //!
    template<typename StateType, typename VectorStateType>
    void get_reactive_scheme( const KineticsConditions<StateType,VectorStateType>& conditions,
                              const VectorStateType& molar_densities,
                              const VectorStateType& h_RT_minus_s_R,
                              VectorStateType& net_rates,
                              VectorStateType& kfwd_const,
                              VectorStateType& kbkwd_const,
                              VectorStateType& kfwd,
                              VectorStateType& kbkwd,
                              VectorStateType& fwd_conc,
                              VectorStateType& bkwd_conc) const;


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

    //! helper function
    //
    // It finds the kinetics model parameter given keywords and reaction:
    //   - nr is the index of the rate to change
    //   - unit is the unit if provided
    //
    // This function is used for both getter and setter.
    void find_kinetics_model_parameter(const unsigned int r,const std::vector<std::string> & keywords, unsigned int & nr, std::string & unit, int & l) const;

    //! helper function
    //
    // It finds the chemical process parameter given keywords and reaction type:
    //   - species is the species whose efficiency we want to change
    //
    // If we want to change anything other than an efficiency, well nothing happens here.
    // This function is used for both getter and setter.
    void find_chemical_process_parameter(ReactionType::Parameters paramChem ,const std::vector<std::string> & keywords, unsigned int & species) const;
     
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
    _reactions.back()->initialize(_reactions.size() - 1);

    return;
  }

  template<typename CoeffType>
  inline
  void ReactionSet<CoeffType>::remove_reaction(unsigned int nr)
  {
     antioch_assert_less(nr,_reactions.size());

     //first clear the memory
     delete _reactions[nr];

     //second, release the spot
     _reactions.erase(_reactions.begin() + nr);
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
  Reaction<CoeffType>& ReactionSet<CoeffType>::reaction(const unsigned int r)
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
      _P0_R(1.0e5/Constants::R_universal<CoeffType>()) //SI
  {
    return;
  }


  template<typename CoeffType>
  inline
  ReactionSet<CoeffType>::~ReactionSet()
  {
    for(unsigned int ir = 0; ir < _reactions.size(); ir++)delete _reactions[ir];
    return;
  }
  
  template<typename CoeffType>
  template<typename StateType, typename VectorStateType, typename VectorReactionsType>
  inline
  void ReactionSet<CoeffType>::compute_reaction_rates ( const KineticsConditions<StateType,VectorStateType>& conditions,
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
    const StateType P0_RT = _P0_R/conditions.T(); // used to transform equilibrium constant from pressure units

    // compute reaction forward rates & other reaction-sized arrays
    for (unsigned int rxn=0; rxn<this->n_reactions(); rxn++)
      {
        net_reaction_rates[rxn] = this->reaction(rxn).compute_rate_of_progress(molar_densities, conditions, P0_RT, h_RT_minus_s_R);
      }
    
    return;
  }

  template<typename CoeffType>
  template<typename StateType, typename VectorStateType, typename VectorReactionsType, typename MatrixReactionsType>
  inline
  void ReactionSet<CoeffType>::compute_reaction_rates_and_derivs( const KineticsConditions<StateType,VectorStateType>& conditions,
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
    const StateType P0_RT = _P0_R/conditions.T(); // used to transform equilibrium constant from pressure units

    // compute reaction forward rates & other reaction-sized arrays
    for (unsigned int rxn=0; rxn<this->n_reactions(); rxn++)
      {
        this->reaction(rxn).compute_rate_of_progress_and_derivatives( molar_densities, _chem_mixture, 
                                                                      conditions, P0_RT, h_RT_minus_s_R, dh_RT_minus_s_R_dT,
                                                                      net_reaction_rates[rxn], 
                                                                      dnet_rate_dT[rxn], 
                                                                      dnet_rate_dX_s[rxn] );
      }
    
    return;
  }


  template<typename CoeffType>
  inline
  unsigned int ReactionSet<CoeffType>::reaction_by_id(const std::string & reaction_id) const
  {
      unsigned int r(0);
      for(r = 0; r < this->n_reactions(); r++)
      {
          if(this->reaction(r).id() == reaction_id)
          break;
      }
      if(r >= this->n_reactions())
      {
        std::string errmsg = "Error: did not find reaction \"" + reaction_id + "\"\nIds are: ";
        for(r = 0; r < this->n_reactions(); r++)
        {
          errmsg += this->reaction(r).id() + ", ";
          if(r%10 == 0)errmsg += "\n"; // a few formatting is nice
        }
        antioch_warning(errmsg);
        antioch_error();
      }

    return r;
  }


  
  template<typename CoeffType>
  inline
  void ReactionSet<CoeffType>::find_kinetics_model_parameter(const unsigned int r,const std::vector<std::string> & keywords, unsigned int & nr, std::string & unit, int & l) const
  {
// now we need to know a few things:
//  if we are in a falloff or duplicate, next keyword is which rate we want, then unit
//  if we are in an elementary photochemistry, next keyword is which index of sigma or lamba we want, then unit
//  else we may have a unit
    switch(this->reaction(r).type())
    {
       case ReactionType::DUPLICATE:
       {
          nr = std::stoi(keywords[1]);           // C++11, throws an exception on error
          if(keywords.size() > 2)unit = keywords[2];
       }
          break;
       case ReactionType::LINDEMANN_FALLOFF:
       case ReactionType::TROE_FALLOFF:
       case ReactionType::LINDEMANN_FALLOFF_THREE_BODY:
       case ReactionType::TROE_FALLOFF_THREE_BODY:
       {
          switch(string_to_kin_enum(keywords[1])) // falloff, if here by error, NOT_FOUND is get (0)
          {
             // This is hard-coded (because we want a vector and not a map)
             // low pressure is the first, high pressure the second
             // see FalloffReaction object
             case KineticsModel::Parameters::LOW_PRESSURE:
             {
                nr = 0;
             }
                break;
             case KineticsModel::Parameters::HIGH_PRESSURE:
             {
                nr = 1;
             }
                break;
              default: // ?
              {
                 antioch_error();
              }
                 break;
           }
           if(keywords.size() > 2)unit = keywords[2]; // unit baby!
       }
          break;
       case ReactionType::ELEMENTARY: // photochem only
       {
          if(this->reaction(r).kinetics_model() == KineticsModel::PHOTOCHEM)
          {
             l = std::stoi(keywords[1]);           // C++11, throws an exception on error
             if(keywords.size() > 2)unit = keywords[2];
          }
       }
          break;
       default:
       {
          if(keywords.size() > 1)unit = keywords[1];  // unit baby!
       }
           break;
    }

    return;
  }

  template<typename CoeffType>
  inline
  void ReactionSet<CoeffType>::find_chemical_process_parameter(ReactionType::Parameters paramChem ,const std::vector<std::string> & keywords, unsigned int & species) const
  {
// there's not really a unit issue here:
//   * efficiencies: unitless
//   * Troe alpha: unitless
//   * Troe T*:   temperature, if it's other than K, you're a bad (like really really bad) scientist, get off my library
//   * Troe T**:  temperature, see above
//   * Troe T***: temperature, do I really need to say anything?

     if(paramChem == ReactionType::Parameters::EFFICIENCIES) // who?
     {
       antioch_assert_greater(keywords.size(),1); // we need a name

       if(!this->chemical_mixture().species_name_map().count(keywords[1]))
                antioch_error(); //who's this?

       species = this->chemical_mixture().species_name_map().at(keywords[1]);
     }

     return;
  }

  template<typename CoeffType>
  inline
  CoeffType ReactionSet<CoeffType>::get_parameter_of_reaction(const std::string & reaction_id, const std::vector<std::string> & keywords) const
  {
     antioch_assert(keywords.size()); // not zero

     CoeffType parameter;
     // when CoeffType is vectorizable, we want
     // to have as little rewriting as possible
     set_zero(parameter);

     unsigned int r = this->reaction_by_id(reaction_id);

     // 1 parse high level
     KineticsModel::Parameters paramKin = string_to_kin_enum(keywords[0]);
     ReactionType::Parameters paramChem = string_to_chem_enum(keywords[0]);

// provide the necessary enum,
// index of reaction rate if kinetics
// index of species if chemical
     if(paramKin != KineticsModel::Parameters::NOT_FOUND)
     {
          // which rate? Duplicate want an unsigned int, falloff a keyword
          unsigned int nr(0); // default is one rate
          std::string unit("SI"); // default internal parameter unit system
          int l(-1); // if vectorized parameter

          this->find_kinetics_model_parameter(r,keywords,nr,unit,l);

          // let's get this
          parameter = reaction(r).get_parameter_of_rate(paramKin, nr, unit);

     }else if(paramChem != ReactionType::Parameters::NOT_FOUND)
     {
          unsigned int species = std::numeric_limits<unsigned int>::max(); // sensible default

          this->find_chemical_process_parameter(paramChem, keywords, species);

          // let's get this
          parameter = reaction(r).get_parameter_of_chemical_process(paramChem, species);
     }else
     {
         antioch_error();
     }

     return parameter;     
 
  }

  template<typename CoeffType>
  template <typename ParamType>
  inline
  void ReactionSet<CoeffType>::set_parameter_of_reaction(const std::string & reaction_id, const std::vector<std::string> & keywords, ParamType value)
  {
     antioch_assert(keywords.size()); // not zero

      // 1 find the reaction
      unsigned int r = this->reaction_by_id(reaction_id);

      // 1 parse high level
      KineticsModel::Parameters paramKin = string_to_kin_enum(keywords[0]);
      ReactionType::Parameters paramChem = string_to_chem_enum(keywords[0]);
// provide the necessary enum,
// index of reaction rate if kinetics
// index of species if chemical
      if(paramKin != KineticsModel::Parameters::NOT_FOUND)
      {
          // which rate? Duplicate want an unsigned int, falloff a keyword
           unsigned int nr(0); // default is one rate
           std::string unit("SI"); // default internal parameter unit system
           int l(-1); // for vectorized parameter

           this->find_kinetics_model_parameter(r,keywords,nr,unit,l);

           // let's set this
           (l < 0)?this->reaction(r).set_parameter_of_rate(paramKin, value, nr, unit):
                   this->reaction(r).set_parameter_of_rate(paramKin, value, nr, l, unit);

      }else if(paramChem != ReactionType::Parameters::NOT_FOUND)
      {
           unsigned int species = std::numeric_limits<unsigned int>::max(); // sensible default

           this->find_chemical_process_parameter(paramChem, keywords, species);

           // let's set this
           this->reaction(r).set_parameter_of_chemical_process(paramChem, value, species);
      }else
      {
          antioch_error();
      }
      
  }

  template<typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void ReactionSet<CoeffType>::print_chemical_scheme( std::ostream& output,
                                                      const KineticsConditions<StateType,VectorStateType>& conditions,
                                                      const VectorStateType& molar_densities,
                                                      const VectorStateType& h_RT_minus_s_R,
                                                      std::vector<VectorStateType>& lossMatrix,
                                                      std::vector<VectorStateType>& prodMatrix,
                                                      std::vector<VectorStateType>& netMatrix ) const
  {

    //filling matrixes
    VectorStateType netRate,kfwd_const,kbkwd_const,kfwd,kbkwd,fwd_conc,bkwd_conc;

    //getting reaction infos
    get_reactive_scheme(conditions,molar_densities,h_RT_minus_s_R,netRate,kfwd_const,kbkwd_const,kfwd,kbkwd,fwd_conc,bkwd_conc);

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
        output << std::setw(15) << kfwd_const[rxn] << "," << std::setw(15) << kbkwd_const[rxn];
      }
    output << std::endl;
    // // conc
    output << std::setw(10) << " ";
    for(unsigned int rxn = 0; rxn < this->n_reactions(); rxn++)
      {
        output << std::setw(15) << fwd_conc[rxn] << "," << std::setw(15) << bkwd_conc[rxn];
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
  void ReactionSet<CoeffType>::get_reactive_scheme( const KineticsConditions<StateType,VectorStateType>& conditions,
                                                    const VectorStateType& molar_densities,
                                                    const VectorStateType& h_RT_minus_s_R,
                                                    VectorStateType& net_rates,
                                                    VectorStateType& kfwd_const,
                                                    VectorStateType& kbkwd_const,
                                                    VectorStateType& kfwd,
                                                    VectorStateType& kbkwd,
                                                    VectorStateType& fwd_conc,
                                                    VectorStateType& bkwd_conc) const
  {
    //!\todo Make these assertions vector-compatible
    // antioch_assert_greater(T, 0.0);

    antioch_assert_equal_to( molar_densities.size(), this->n_species() );
    antioch_assert_equal_to( h_RT_minus_s_R.size(), this->n_species() );

    // useful constants
    const StateType P0_RT = _P0_R/conditions.T(); // used to transform equilibrium constant from pressure units

    net_rates.resize(this->n_reactions(),0);
    kfwd_const.resize(this->n_reactions(),0);
    kfwd.resize(this->n_reactions(),0);
    kbkwd_const.resize(this->n_reactions(),0);
    kbkwd.resize(this->n_reactions(),0);
    fwd_conc.resize(this->n_reactions(),1);
    bkwd_conc.resize(this->n_reactions(),1);

    // compute reaction forward rates & other reaction-sized arrays
    for (unsigned int rxn=0; rxn<this->n_reactions(); rxn++)
      {
        const Reaction<CoeffType>& reaction = this->reaction(rxn);
        kfwd_const[rxn] = reaction.compute_forward_rate_coefficient(molar_densities,conditions);
        kfwd[rxn] = kfwd_const[rxn];

        for (unsigned int r=0; r<reaction.n_reactants(); r++)
          {
            fwd_conc[rxn] *= pow( molar_densities[reaction.reactant_id(r)],
                              reaction.reactant_partial_order(r));
          }
        kfwd[rxn] *= fwd_conc[rxn];

        if(reaction.reversible())
        {
          StateType keq = reaction.equilibrium_constant( P0_RT, h_RT_minus_s_R );

          kbkwd_const[rxn] = kfwd_const[rxn]/keq;
          kbkwd[rxn] = kbkwd_const[rxn];
          for (unsigned int p=0; p<reaction.n_products(); p++)
            {
              bkwd_conc[rxn] *= pow( molar_densities[reaction.product_id(p)],
                                  reaction.product_partial_order(p));
            }
          kbkwd[rxn] *= bkwd_conc[rxn];
        }

        net_rates[rxn] = kfwd[rxn] - kbkwd[rxn];

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
           << (this->reaction(r)) << "\n";
      }

    return;
  }
 
} // end namespace Antioch

#endif // ANTIOCH_REACTION_SET_H
