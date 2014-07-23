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

#ifndef ANTIOCH_REACTION_H
#define ANTIOCH_REACTION_H

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/cmath_shims.h"
#include "antioch/kinetics_type.h"
#include "antioch/constant_rate.h"
#include "antioch/hercourtessen_rate.h"
#include "antioch/berthelot_rate.h"
#include "antioch/arrhenius_rate.h"
#include "antioch/berthelothercourtessen_rate.h"
#include "antioch/kooij_rate.h"
#include "antioch/vanthoff_rate.h"
#include "antioch/photochemical_rate.h"
#include "antioch/reaction_enum.h"
#include "antioch/chemical_mixture.h"
#include "antioch/kinetics_conditions.h"

//C++
#include <string>
#include <vector>
#include <iostream>

namespace Antioch
{
  // Forward declarations
  template <typename CoeffType>
  class ElementaryReaction;

  template <typename CoeffType>
  class DuplicateReaction;

  template <typename CoeffType>
  class ThreeBodyReaction;

  template <typename CoeffType,typename FalloffType>
  class FalloffReaction;

  template <typename CoeffType>
  class LindemannFalloff;

  template <typename CoeffType>
  class TroeFalloff;

  //!A single reaction mechanism. 
  /*!\class Reaction
   *
   * A reaction is characterized by the rate constant \f$k(T,[M])\f$, which can
   * be decomposed into a chemical process and a kinetics model:
   * \f[
   *  k(T,[M]) = \chi([M],\alpha(T))
   * \f]
   * with \f$\chi\f$ the chemical process and \f$\alpha\f$ the kinetics model.
   *
   This virtual base is derived for the following processes:
   - elementary process (ElementaryReaction),
   - duplicate process (DuplicateReaction),
   - three body process (ThreeBodyReaction)
   - falloff processes (FalloffReaction) with:
      - Lindemann falloff (LindemannFalloff),
      - Troe falloff (TroeFalloff).

   This class encapsulates a kinetics model.  The choosable kinetics models are
   - Constant (ConstantRate),
   - Hercourt Hessen (HercourtEssenRate),
   - Berthelot  (BerthelotRate),
   - Arrhenius  (ArrheniusRate),
   - Berthelot Hercourt Hessen  (BerthelotHercourtEssenRate),
   - Kooij, or modified Arrhenius  (KooijRate),
   - Van't Hoff  (VantHoffRate).

   By default, we choose a reversible ElementaryReaction with a KooijRate kinetics model.
  */
  template<typename CoeffType=double, typename VectorCoeffType = std::vector<CoeffType> >
  class Reaction
  {
  public:

    //! Construct a single reaction mechanism.
    Reaction( const unsigned int n_species, const std::string &equation,
              const bool &reversible = true,
              const ReactionType::ReactionType type = ReactionType::ELEMENTARY,
              const KineticsModel::KineticsModel kin = KineticsModel::KOOIJ);
    
    virtual ~Reaction();

    unsigned int n_species() const;
    
    //! \returns the equation for this reaction.
    std::string equation() const;

    /*! Type of reaction.
     *  reversible reactions are considered.
     */
    ReactionType::ReactionType type() const;

    /*! Set the type of reaction.
     * reversible reactions are considered.
     */
    void set_type( const ReactionType::ReactionType type);

    /*! Set the reversibility of reaction.
     */
    void set_reversibility( const bool reversible);

    //! Model of kinetics.
    KineticsModel::KineticsModel kinetics_model() const;

    //! Set the model of kinetics.
    void set_kinetics_model( const KineticsModel::KineticsModel kin);
    
    bool initialized() const;


    /*! \return the reversibility state of reaction.
     */
    bool reversible() const;

    //! \returns the number of reactants.
    unsigned int n_reactants() const;

    //! \returns the number of products.
    unsigned int n_products() const;

    //! \returns the name of the \p r th reactant.
    const std::string& reactant_name(const unsigned int r) const;

    //! \returns the name of the \p p th product.
    const std::string& product_name(const unsigned int p) const;

    //!
    unsigned int reactant_id(const unsigned int r) const;

    //!
    unsigned int product_id(const unsigned int p) const;

    //!
    unsigned int reactant_stoichiometric_coefficient(const unsigned int r) const;

    //!
    unsigned int product_stoichiometric_coefficient(const unsigned int p) const;

    //!
    void add_reactant( const std::string &name,
                       const unsigned int r_id,
                       const unsigned int stoichiometric_coeff);

    //!
    void add_product( const std::string &name,
                      const unsigned int p_id,
                      const unsigned int stoichiometric_coeff);

    //!
    void set_efficiency( const std::string &,
                         const unsigned int s,
                         const CoeffType efficiency);

    //!
    CoeffType efficiency( const unsigned int s) const;

    //! Computes derived quantities.
    void initialize(unsigned int index = 0);

    //!
    int gamma() const;
    
    //!
    template <typename StateType, typename VectorStateType>
    StateType equilibrium_constant( const StateType& P0_RT,
                                    const VectorStateType& h_RT_minus_s_R ) const;

    //!
    template <typename StateType, typename VectorStateType>
    void equilibrium_constant_and_derivative( const StateType& T,
                                              const StateType& P0_RT,
                                              const VectorStateType& h_RT_minus_s_R,
                                              const VectorStateType& ddT_h_RT_minus_s_R,
                                              StateType& keq,
                                              StateType& dkeq_dT) const;

    //!
    template <typename StateType, typename VectorStateType>
    StateType compute_forward_rate_coefficient( const VectorStateType& molar_densities,
                                                const KineticsConditions<StateType,VectorStateType>& conditions) const;
    
    //!
    template <typename StateType, typename VectorStateType>
    void compute_forward_rate_coefficient_and_derivatives( const VectorStateType& molar_densities,
                                                           const KineticsConditions<StateType,VectorStateType>& conditions, 
                                                           StateType& kfwd, 
                                                           StateType& dkfwd_dT,
                                                           VectorStateType& dkfwd_dX) const;
    ////
    template <typename StateType, typename VectorStateType>
    StateType compute_rate_of_progress( const VectorStateType& molar_densities,
                                        const KineticsConditions<StateType,VectorStateType>& conditions,  
                                        const StateType& P0_RT,  
                                        const VectorStateType& h_RT_minus_s_R) const;

    template <typename StateType, typename VectorStateType>
    void compute_rate_of_progress_and_derivatives( const VectorStateType &molar_densities,
                                                   const ChemicalMixture<CoeffType>& chem_mixture,
                                                   const KineticsConditions<StateType,VectorStateType>& conditions,
                                                   const StateType &P0_RT,
                                                   const VectorStateType &h_RT_minus_s_R,
                                                   const VectorStateType &dh_RT_minus_s_R_dT,
                                                   StateType& net_reaction_rate,
                                                   StateType& dnet_rate_dT,
                                                   VectorStateType& dnet_rate_dX_s ) const;

    //! Return const reference to the forward rate object
    const KineticsType<CoeffType,VectorCoeffType>& forward_rate(unsigned int ir = 0) const;

    //! Return writeable reference to the forward rate object
    KineticsType<CoeffType,VectorCoeffType>& forward_rate(unsigned int ir = 0);

    //! Add a forward rate object
    void add_forward_rate(KineticsType<CoeffType,VectorCoeffType> *rate);

    //! Swap two forward rates object
    void swap_forward_rates(unsigned int irate, unsigned int jrate);

    //! Return the number of rate constant objects
    unsigned int n_rate_constants() const;

    //! Formatted print, by default to \p std::cout.
    void print(std::ostream& os = std::cout) const;

    //! Formatted print.
    friend std::ostream& operator << (std::ostream& os, const Reaction &rxn)
    {
      rxn.print(os);
      return os;
    }

  protected:

    unsigned int _n_species;
    std::string _equation;
    std::vector<std::string> _reactant_names;
    std::vector<std::string> _product_names;
    std::vector<unsigned int> _reactant_ids;
    std::vector<unsigned int> _product_ids;
    std::vector<unsigned int> _reactant_stoichiometry;
    std::vector<unsigned int> _product_stoichiometry;
    std::vector<CoeffType> _efficiencies;
    std::vector<unsigned int> _species_reactant_stoichiometry;
    std::vector<unsigned int> _species_product_stoichiometry;
    std::vector<int>          _species_delta_stoichiometry;
    int _gamma;
    bool _initialized;
    bool _reversible;
    ReactionType::ReactionType _type;
    KineticsModel::KineticsModel _kintype;

    //! The forward reaction rate modified Arrhenius form.
    std::vector<KineticsType<CoeffType,VectorCoeffType>* > _forward_rate;

  private:
    Reaction();

  };

  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType, typename VectorCoeffType>
  inline
  unsigned int Reaction<CoeffType,VectorCoeffType>::n_species() const
  {
    return _n_species;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  std::string Reaction<CoeffType,VectorCoeffType>::equation() const
  {
    return _equation;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  ReactionType::ReactionType Reaction<CoeffType,VectorCoeffType>::type() const
  {
    return _type;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  void Reaction<CoeffType,VectorCoeffType>::set_type( const ReactionType::ReactionType type)
  {
    _type = type;
    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  void Reaction<CoeffType,VectorCoeffType>::set_reversibility( const bool reversible)
  {
    _reversible = reversible;
    return;
  }

  template<typename CoeffType,typename VectorCoeffType>
  inline
  KineticsModel::KineticsModel Reaction<CoeffType,VectorCoeffType>::kinetics_model() const
  {
    return _kintype;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  void Reaction<CoeffType,VectorCoeffType>::set_kinetics_model( const KineticsModel::KineticsModel kin)
  {
    _kintype = kin;
    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  bool Reaction<CoeffType,VectorCoeffType>::initialized() const
  {
    return _initialized;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  bool Reaction<CoeffType,VectorCoeffType>::reversible() const
  {
    return _reversible;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  void Reaction<CoeffType,VectorCoeffType>::add_forward_rate(KineticsType<CoeffType,VectorCoeffType> *rate)
  {
    _forward_rate.push_back(rate);
    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  unsigned int Reaction<CoeffType,VectorCoeffType>::n_rate_constants() const
  {
    return _forward_rate.size();
  }
    
  template<typename CoeffType, typename VectorCoeffType>
  inline
  void Reaction<CoeffType,VectorCoeffType>::swap_forward_rates(unsigned int irate, unsigned int jrate)
  {
    antioch_assert_less(irate,_forward_rate.size());
    antioch_assert_less(jrate,_forward_rate.size());
    
    KineticsType<CoeffType>* rate_tmp = _forward_rate[jrate];
    _forward_rate[jrate] = _forward_rate[irate];
    _forward_rate[irate] = rate_tmp;

    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  unsigned int Reaction<CoeffType,VectorCoeffType>::n_reactants () const
  {
    antioch_assert_less(_reactant_ids.size(), this->n_species());
    antioch_assert_equal_to(_reactant_ids.size(), _reactant_stoichiometry.size());
    antioch_assert_equal_to(_reactant_ids.size(), _reactant_names.size());
    return _reactant_ids.size();
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  unsigned int Reaction<CoeffType,VectorCoeffType>::n_products() const
  {
    antioch_assert_less(_product_ids.size(), this->n_species());
    antioch_assert_equal_to(_product_ids.size(), _product_stoichiometry.size());
    antioch_assert_equal_to(_product_ids.size(),  _product_names.size());
    return _product_ids.size();
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  const std::string& Reaction<CoeffType,VectorCoeffType>::reactant_name(const unsigned int r) const
  {       
    antioch_assert_less(r, _reactant_names.size()); 
    return _reactant_names[r]; 
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  const std::string& Reaction<CoeffType,VectorCoeffType>::product_name(const unsigned int p) const
  { 
    antioch_assert_less(p, _product_names.size()); 
    return _product_names[p]; 
  }
  
  template<typename CoeffType, typename VectorCoeffType>
  inline
  unsigned int Reaction<CoeffType,VectorCoeffType>::reactant_id(const unsigned int r) const
  { 
    antioch_assert_less(r, _reactant_ids.size()); 
    antioch_assert_less(_reactant_ids[r], this->n_species());
    return _reactant_ids[r]; 
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  unsigned int Reaction<CoeffType,VectorCoeffType>::product_id(const unsigned int p) const
  { 
    antioch_assert_less(p, _product_ids.size()); 
    antioch_assert_less(_product_ids[p], this->n_species());
    return _product_ids[p]; 
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  unsigned int Reaction<CoeffType,VectorCoeffType>::reactant_stoichiometric_coefficient(const unsigned int r) const
  {      
    antioch_assert_less(r, _reactant_stoichiometry.size());
    antioch_assert_less(_reactant_ids[r], this->n_species());
    return _reactant_stoichiometry[r];
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  unsigned int Reaction<CoeffType,VectorCoeffType>::product_stoichiometric_coefficient(const unsigned int p) const
  {
    antioch_assert_less(p, _product_stoichiometry.size());
    antioch_assert_less(_product_ids[p], this->n_species());
    return _product_stoichiometry[p];
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  void Reaction<CoeffType,VectorCoeffType>::add_reactant (const std::string &name,
                                          const unsigned int r_id,
                                          const unsigned int stoichiometric_coeff)
  {
    antioch_assert_less(r_id, this->n_species());
    _reactant_names.push_back(name);
    _reactant_ids.push_back(r_id);
    _reactant_stoichiometry.push_back(stoichiometric_coeff);
    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  void Reaction<CoeffType,VectorCoeffType>::add_product (const std::string &name,
                                         const unsigned int p_id,
                                         const unsigned int stoichiometric_coeff)
  {
    antioch_assert_less(p_id, this->n_species());
    _product_names.push_back(name);
    _product_ids.push_back(p_id);
    _product_stoichiometry.push_back(stoichiometric_coeff);
    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  void Reaction<CoeffType,VectorCoeffType>::set_efficiency (const std::string &,
                                            const unsigned int s,
                                            const CoeffType efficiency)
  {
    antioch_assert_less(s, this->n_species());
    antioch_assert_less(s, _efficiencies.size());
    antioch_assert_equal_to(_type, ReactionType::THREE_BODY);
    _efficiencies[s] = efficiency;
    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  CoeffType Reaction<CoeffType,VectorCoeffType>::efficiency( const unsigned int s ) const
  {
    antioch_assert_less(s, _efficiencies.size());
    antioch_assert_equal_to(_type, ReactionType::THREE_BODY);
    return _efficiencies[s];
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  int Reaction<CoeffType,VectorCoeffType>::gamma () const
  {
    return _gamma;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  const KineticsType<CoeffType,VectorCoeffType>& Reaction<CoeffType,VectorCoeffType>::forward_rate(unsigned int ir) const
  {
    antioch_assert_less(ir,_forward_rate.size());
    return *_forward_rate[ir];
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  KineticsType<CoeffType,VectorCoeffType>& Reaction<CoeffType,VectorCoeffType>::forward_rate(unsigned int ir)
  {
    antioch_assert_less(ir,_forward_rate.size());
    return *_forward_rate[ir];
  }


  template<typename CoeffType, typename VectorCoeffType>
  inline
  Reaction<CoeffType,VectorCoeffType>::Reaction( const unsigned int n_species,
                                 const std::string &equation,
                                 const bool &reversible,
                                 const ReactionType::ReactionType type,
                                 const KineticsModel::KineticsModel kin)
    : _n_species(n_species),
      _equation(equation),
      _gamma(0),
      _initialized(false),
      _reversible(reversible),
      _type(type),
      _kintype(kin)
  {
    _efficiencies.resize(_n_species); 
    std::fill (_efficiencies.begin(), _efficiencies.end(), 1.);
  }


  template<typename CoeffType, typename VectorCoeffType>
  inline
  Reaction<CoeffType,VectorCoeffType>::~Reaction()
  {
    for(unsigned int ir = 0; ir < _forward_rate.size(); ir++)
      {
        delete _forward_rate[ir];
      }

    return;
  }


  template<typename CoeffType, typename VectorCoeffType>
  inline
  void Reaction<CoeffType,VectorCoeffType>::initialize(unsigned int index)
  {
    // Stoichiometric coefficients, by species id
    {
      _species_reactant_stoichiometry.resize(this->n_species());
      _species_product_stoichiometry.resize(this->n_species());
      _species_delta_stoichiometry.resize(this->n_species());

      std::fill( _species_reactant_stoichiometry.begin(),
                 _species_reactant_stoichiometry.end(),
                 0 );

      std::fill( _species_product_stoichiometry.begin(),
                 _species_product_stoichiometry.end(),
                 0);
    }
    
    for (unsigned int r=0; r< this->n_reactants(); r++)
      {
        _species_reactant_stoichiometry[this->reactant_id(r)] =
          this->reactant_stoichiometric_coefficient(r);
      }
    
    for (unsigned int p=0; p < this->n_products(); p++)
      {
        _species_product_stoichiometry[this->product_id(p)] =
          this->product_stoichiometric_coefficient(p);
      }
    
    // find the delta stoichiometric coefficient for each species,
    // and the sum of the deltas 
    _gamma = 0;
    for (unsigned int s=0; s<this->n_species(); s++)
      {
        _species_delta_stoichiometry[s] =
          ( _species_product_stoichiometry[s] -
            _species_reactant_stoichiometry[s] );
        
        _gamma += _species_delta_stoichiometry[s];
      }

     // gives kinetics object index in reaction set
     for(typename std::vector<KineticsType<CoeffType,VectorCoeffType>* >::iterator it = _forward_rate.begin();
                it != _forward_rate.end(); it++)
     {
        (*it)->set_index(index);
     }
   
    // set initialization flag
    _initialized = true;
  }


  template<typename CoeffType, typename VectorCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  StateType Reaction<CoeffType,VectorCoeffType>::equilibrium_constant( const StateType& P0_RT,
                                                       const VectorStateType& h_RT_minus_s_R ) const
  {
    using std::exp;

    antioch_assert( this->initialized() );
    //!\todo Make this assertion vector-compatible
    // antioch_assert_greater( P0_RT, 0.0 );
    antioch_assert_greater( h_RT_minus_s_R.size(), 0 );
    antioch_assert_equal_to( h_RT_minus_s_R.size(), this->n_species() );
    antioch_assert_equal_to( _species_delta_stoichiometry.size(), this->n_species() );

    StateType exppower = -( static_cast<CoeffType>(_species_delta_stoichiometry[0])*
                            h_RT_minus_s_R[0] );

    for (unsigned int s=1; s < this->n_species(); s++)
      {
        exppower += -( static_cast<CoeffType>(_species_delta_stoichiometry[s])*
                       h_RT_minus_s_R[s] );
      }

    return ant_pow( P0_RT, static_cast<CoeffType>(this->gamma()) )*exp(exppower);
  }


  template<typename CoeffType, typename VectorCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void Reaction<CoeffType,VectorCoeffType>::equilibrium_constant_and_derivative( const StateType& T,
                                                                 const StateType& P0_RT,
                                                                 const VectorStateType& h_RT_minus_s_R,
                                                                 const VectorStateType& ddT_h_RT_minus_s_R,
                                                                 StateType& keq,
                                                                 StateType& dkeq_dT) const
  {
    antioch_assert(this->initialized());

    //!\todo Make these assertions vector-compatible
    // antioch_assert_greater( P0_RT, 0.0 );
    // antioch_assert_greater( T, 0.0 );
    antioch_assert_greater( h_RT_minus_s_R.size(), 0 );
    antioch_assert_equal_to( h_RT_minus_s_R.size(), this->n_species() );
    antioch_assert_equal_to( ddT_h_RT_minus_s_R.size(), this->n_species() );
    antioch_assert_equal_to( _species_delta_stoichiometry.size(), this->n_species() );

    // get the equilibrium constant
    keq = this->equilibrium_constant( P0_RT, h_RT_minus_s_R );

    StateType ddT_exppower = -( static_cast<CoeffType>(_species_delta_stoichiometry[0])*
                                ddT_h_RT_minus_s_R[0] );

    for (unsigned int s=1; s<this->n_species(); s++)
      ddT_exppower += -( static_cast<CoeffType>(_species_delta_stoichiometry[s])*
                         ddT_h_RT_minus_s_R[s] );

    // compute its derivative
    dkeq_dT = keq*(-static_cast<CoeffType>(this->gamma())/T + ddT_exppower);

    return;
  }


  template<typename CoeffType, typename VectorCoeffType>
  inline
  void Reaction<CoeffType,VectorCoeffType>::print( std::ostream& os ) const
  {
    os << "# Gas-Phase Reaction \"" << _equation << "\":\n";
    if (this->n_species())
      {
        os << "#   reactants: ";
        for (unsigned int r=0; r<this->n_reactants(); r++)
          os << this->reactant_name(r) << ":"
             << this->reactant_stoichiometric_coefficient(r) << " ";
        os << "\n"
           << "#   products:  ";
        for (unsigned int p=0; p<this->n_products(); p++)
          os << this->product_name(p) << ":"
             << this->product_stoichiometric_coefficient(p) << " ";
      }
    os << "\n#   Chemical process: " << _type;
    os << "\n#   Kinetics model: "   << _kintype;
    for(unsigned int ir = 0; ir < _forward_rate.size(); ir++)
      {
        os << "\n#   forward rate eqn: " << *_forward_rate[ir];
      }

    if (_type == ReactionType::THREE_BODY)
      {
        os << "\n#   efficiencies: ";
        for (unsigned int s=0; s<this->n_species(); s++)
          os << s << ":" << this->efficiency(s) << " ";
      }
    os << "\n#";
    return;
  }


  template<typename CoeffType, typename VectorCoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  StateType Reaction<CoeffType,VectorCoeffType>::compute_forward_rate_coefficient( const VectorStateType& molar_densities,
                                                                   const KineticsConditions<StateType,VectorStateType>& conditions) const
  {
    switch(_type)
      {
      case(ReactionType::ELEMENTARY):
        {
          return (static_cast<const ElementaryReaction<CoeffType>*>(this))->compute_forward_rate_coefficient(molar_densities,conditions);
        }
        break;

      case(ReactionType::DUPLICATE):
        {
          return (static_cast<const DuplicateReaction<CoeffType>*>(this))->compute_forward_rate_coefficient(molar_densities,conditions);
        }
        break;

      case(ReactionType::THREE_BODY):
        {
          return (static_cast<const ThreeBodyReaction<CoeffType>*>(this))->compute_forward_rate_coefficient(molar_densities,conditions);
        }
        break;

      case(ReactionType::LINDEMANN_FALLOFF):
        {
          return (static_cast<const FalloffReaction<CoeffType,LindemannFalloff<CoeffType> >*>(this))->compute_forward_rate_coefficient(molar_densities,conditions);
        }
        break;

      case(ReactionType::TROE_FALLOFF):
        {
          return (static_cast<const FalloffReaction<CoeffType,TroeFalloff<CoeffType> >*>(this))->compute_forward_rate_coefficient(molar_densities,conditions);
        }
        break;

      default:
        {
          antioch_error();
        }
      } // switch(_type)

    // Dummy
    return zero_clone(conditions.T());
  }

  template<typename CoeffType, typename VectorCoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  void Reaction<CoeffType,VectorCoeffType>::compute_forward_rate_coefficient_and_derivatives( const VectorStateType& molar_densities,
                                                                              const KineticsConditions<StateType,VectorStateType>& conditions, 
                                                                              StateType& kfwd, 
                                                                              StateType& dkfwd_dT,
                                                                              VectorStateType& dkfwd_dX) const
  {
    switch(_type)
      {
      case(ReactionType::ELEMENTARY):
        {
          (static_cast<const ElementaryReaction<CoeffType>*>(this))->compute_forward_rate_coefficient_and_derivatives(molar_densities,conditions,kfwd,dkfwd_dT,dkfwd_dX);
        }
        break;

      case(ReactionType::DUPLICATE):
        {
          (static_cast<const DuplicateReaction<CoeffType>*>(this))->compute_forward_rate_coefficient_and_derivatives(molar_densities,conditions,kfwd,dkfwd_dT,dkfwd_dX);
        }
        break;

      case(ReactionType::THREE_BODY):
        {
          (static_cast<const ThreeBodyReaction<CoeffType>*>(this))->compute_forward_rate_coefficient_and_derivatives(molar_densities,conditions,kfwd,dkfwd_dT,dkfwd_dX);
        }
        break;

      case(ReactionType::LINDEMANN_FALLOFF):
        {
          (static_cast<const FalloffReaction<CoeffType,LindemannFalloff<CoeffType> >*>(this))->compute_forward_rate_coefficient_and_derivatives(molar_densities,conditions,kfwd,dkfwd_dT,dkfwd_dX);
        }
        break;

      case(ReactionType::TROE_FALLOFF):
        {
          (static_cast<const FalloffReaction<CoeffType,TroeFalloff<CoeffType> >*>(this))->compute_forward_rate_coefficient_and_derivatives(molar_densities,conditions,kfwd,dkfwd_dT,dkfwd_dX);
        }
        break;

      default:
        {
          antioch_error();
        }

      } // switch(type)

    return;
  }

  //kfwd *prod_r [R]^nu_r - kbkwd * prod_p [P]^nu_p ( = - 1/nu_r d[R]/dt)
  template<typename CoeffType, typename VectorCoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  StateType Reaction<CoeffType,VectorCoeffType>::compute_rate_of_progress( const VectorStateType& molar_densities,
                                                           const KineticsConditions<StateType,VectorStateType>& conditions,  
                                                           const StateType& P0_RT,  
                                                           const VectorStateType& h_RT_minus_s_R) const
  {
    StateType kfwd = this->compute_forward_rate_coefficient(molar_densities,conditions);
    StateType kfwd_times_reactants = kfwd;

    // Rfwd
    for (unsigned int ro=0; ro < this->n_reactants(); ro++)
      {
        kfwd_times_reactants     *= 
          ant_pow( molar_densities[this->reactant_id(ro)],
            static_cast<int>(this->reactant_stoichiometric_coefficient(ro)) );
      }

    StateType kbkwd_times_products = Antioch::zero_clone(kfwd_times_reactants);
    if(_reversible)
    {

      StateType Keq = this->equilibrium_constant( P0_RT, h_RT_minus_s_R );
      kbkwd_times_products = kfwd/Keq;

      // Rbkwd 
      for (unsigned int po=0; po< this->n_products(); po++)
        {
          kbkwd_times_products     *= 
            ant_pow( molar_densities[this->product_id(po)],
              static_cast<int>(this->product_stoichiometric_coefficient(po)) );
         
        }
    }

    return kfwd_times_reactants - kbkwd_times_products;

  }

  template<typename CoeffType, typename VectorCoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  void Reaction<CoeffType,VectorCoeffType>::compute_rate_of_progress_and_derivatives( const VectorStateType &molar_densities,
                                                                      const ChemicalMixture<CoeffType>& chem_mixture,
                                                                      const KineticsConditions<StateType,VectorStateType>& conditions,
                                                                      const StateType &P0_RT,
                                                                      const VectorStateType &h_RT_minus_s_R,
                                                                      const VectorStateType &dh_RT_minus_s_R_dT,
                                                                      StateType& net_reaction_rate,
                                                                      StateType& dnet_rate_dT,
                                                                      VectorStateType& dnet_rate_dX_s ) const
  {
    antioch_assert_equal_to (molar_densities.size(), this->n_species());

// First the forward component, if reversible, compute and add the backward component

    StateType kfwd = Antioch::zero_clone(conditions.T());
    StateType dkfwd_dT = Antioch::zero_clone(conditions.T());
    VectorStateType dkfwd_dX_s = Antioch::zero_clone(molar_densities);
    VectorStateType dkbkwd_dX_s = Antioch::zero_clone(molar_densities);

    this->compute_forward_rate_coefficient_and_derivatives(molar_densities, conditions, kfwd, dkfwd_dT ,dkfwd_dX_s);

    // If users want to use valarrays, then the output reference sizes
    // had better already match the input value sizes...
    StateType dRfwd_dT = Antioch::zero_clone(kfwd);

    // We need to construct using an input StateType argument if we
    // want StateType==valarray to have the right sizes
    // valarray compatibility makes this a bit redundant, but not much
    // worse than the previous version
    /*! \todo Should we make this work arrays that get passed in so we aren't allocating/deallocating here? */
    VectorStateType dRfwd_dX_s(this->n_species(), kfwd);

    Antioch::set_zero(dRfwd_dX_s);

    // pre-fill the participating species partials with the rates
    for (unsigned int r=0; r< this->n_reactants(); r++)
      {
        dRfwd_dX_s[this->reactant_id(r)] = kfwd;
      }
    
    //init
    StateType facfwd = constant_clone(conditions.T(),1);
    dRfwd_dT = dkfwd_dT;

    // Rfwd & derivatives
    for (unsigned int ro=0; ro < this->n_reactants(); ro++)
      {
        const StateType val = 
          ant_pow( molar_densities[this->reactant_id(ro)],
            static_cast<int>(this->reactant_stoichiometric_coefficient(ro)) );
          
        const StateType dval = 
          ( static_cast<CoeffType>(this->reactant_stoichiometric_coefficient(ro))*
            ant_pow( molar_densities[this->reactant_id(ro)],
              static_cast<int>(this->reactant_stoichiometric_coefficient(ro))-1 ) 
            );

        facfwd   *= val;
        dRfwd_dT *= val;

        for (unsigned int ri=0; ri<this->n_reactants(); ri++)
          {
            dRfwd_dX_s[this->reactant_id(ri)] *= (ri == ro) ? dval : val;
          }
      }

    for (unsigned int s = 0; s < this->n_species(); s++)
      {
        dRfwd_dX_s[s] += facfwd * dkfwd_dX_s[s];
      }
        
    net_reaction_rate = facfwd * kfwd;

    dnet_rate_dT = dRfwd_dT;

    for (unsigned int s = 0; s < this->n_species(); s++)
      {
        dnet_rate_dX_s[s] = dRfwd_dX_s[s];
      }

    if(_reversible)
    {

    //backward to be computed and added

      StateType keq = Antioch::zero_clone(conditions.T());
      StateType dkeq_dT = Antioch::zero_clone(conditions.T());

      equilibrium_constant_and_derivative( conditions.T(), P0_RT, h_RT_minus_s_R,
                                           dh_RT_minus_s_R_dT,
                                           keq, dkeq_dT );

      const StateType kbkwd = kfwd/keq;
      const StateType dkbkwd_dT = (dkfwd_dT - kbkwd*dkeq_dT)/keq;
      for(unsigned int s = 0; s < this->n_species(); s++)
        {
          dkbkwd_dX_s[s] = dkfwd_dX_s[s]/keq;
        }

      // If users want to use valarrays, then the output reference sizes
      // had better already match the input value sizes...
      StateType dRbkwd_dT = Antioch::zero_clone(dkbkwd_dT);

      // We need to construct using an input StateType argument if we
      // want StateType==valarray to have the right sizes
      // valarray compatibility makes this a bit redundant, but not much
      // worse than the previous version
      /*! \todo Should we make this work arrays that get passed in so we aren't allocating/deallocating here? */
      VectorStateType dRbkwd_dX_s(this->n_species(), kbkwd);

      Antioch::set_zero(dRbkwd_dX_s);

      // pre-fill the participating species partials with the rates
      for (unsigned int p=0; p < this->n_products(); p++)
        {
          dRbkwd_dX_s[this->product_id(p)] = kbkwd;
        }

      //init
      StateType facbkwd = constant_clone(conditions.T(),1);
      dRbkwd_dT = dkbkwd_dT;

      // Rbkwd & derivatives
      for (unsigned int po=0; po< this->n_products(); po++)
        {
          const StateType val = 
            ant_pow( molar_densities[this->product_id(po)],
              static_cast<int>(this->product_stoichiometric_coefficient(po)) );
          
          const StateType dval = 
            ( static_cast<CoeffType>(this->product_stoichiometric_coefficient(po))*
              ant_pow( molar_densities[this->product_id(po)],
                static_cast<int>(this->product_stoichiometric_coefficient(po))-1 )
              );
        
          facbkwd   *= val;
          dRbkwd_dT *= val;

          for (unsigned int pi=0; pi<this->n_products(); pi++)
            {
              dRbkwd_dX_s[this->product_id(pi)] *= (pi == po) ? dval : val;
            }
        
        }

      for (unsigned int s = 0; s < this->n_species(); s++)
        {
          dRbkwd_dX_s[s] += facbkwd * dkbkwd_dX_s[s];
        }

      net_reaction_rate -= facbkwd * kbkwd;

      dnet_rate_dT -= dRbkwd_dT;

      for (unsigned int s = 0; s < this->n_species(); s++)
        {
          dnet_rate_dX_s[s] -= dRbkwd_dX_s[s];
        }

    } //end of if(_reversible) condition

    return;
  }
    
} // namespace Antioch

#endif // ANTIOCH_REACTION_H
