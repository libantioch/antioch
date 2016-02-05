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
#include "antioch/kinetics_parsing.h" // reset_parameter_of_rate

//C++
#include <string>
#include <vector>
#include <iostream>
#include <limits>

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

  template <typename CoeffType,typename FalloffType>
  class FalloffThreeBodyReaction;

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
   - falloff processes (FalloffReaction) and falloff three-body processes (FalloffThreeBodyReaction) with:
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
    const std::string & equation() const;

    //! \returns the reaction id.
    const std::string & id() const;

    //! set the reaction id.
    void set_id(const std::string & id);

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

    //! reset a parameter from the rate constant
    void set_parameter_of_rate(KineticsModel::Parameters parameter, CoeffType new_value, unsigned int n_reaction = 0, const std::string & unit = "SI");

    //! reset a parameter from the rate constant, vector parameters
    void set_parameter_of_rate(KineticsModel::Parameters parameter, CoeffType new_value, unsigned int n_reaction, int l, const std::string & unit = "SI");


    //! get a parameter from the rate constant
    CoeffType get_parameter_of_rate(KineticsModel::Parameters parameter, unsigned int n_reaction = 0, const std::string & unit = "SI") const;

    //! get a parameter from the rate constant, vectorized version
    CoeffType get_parameter_of_rate(KineticsModel::Parameters parameter, unsigned int n_reaction, const std::string &  unit, int l) const;

    //! reset a parameter from the chemical process
    void set_parameter_of_chemical_process(ReactionType::Parameters parameter, CoeffType new_value, unsigned int species = std::numeric_limits<unsigned int>::max() );

    //! get a parameter from the chemical process
    CoeffType get_parameter_of_chemical_process(ReactionType::Parameters parameter, unsigned int species = std::numeric_limits<unsigned int>::max() ) const;

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
    CoeffType reactant_partial_order(const unsigned int r) const;

    //!
    CoeffType product_partial_order(const unsigned int p) const;

    //!
    void add_reactant( const std::string &name,
                       const unsigned int r_id,
                       const unsigned int stoichiometric_coeff,
                       const CoeffType partial_order = std::numeric_limits<CoeffType>::infinity());// what test could be reliable?

    //!
    void add_product( const std::string &name,
                      const unsigned int p_id,
                      const unsigned int stoichiometric_coeff,
                      const CoeffType partial_order = std::numeric_limits<CoeffType>::infinity()); // what test could be reliable?

    //!
    void clear_reactant();

    //!
    void clear_product();

    //!
    void set_efficiency( const std::string &,
                         const unsigned int s,
                         const CoeffType efficiency);

    //!
    CoeffType get_efficiency( const unsigned int s) const;

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

    // Deprecated API for backwards compatibility
    template <typename StateType, typename VectorStateType>
    StateType compute_forward_rate_coefficient( const VectorStateType& molar_densities,
                                                const StateType& temp)
    const;

    //!
    //!
    template <typename StateType, typename VectorStateType>
    void compute_forward_rate_coefficient_and_derivatives( const VectorStateType& molar_densities,
                                                           const KineticsConditions<StateType,VectorStateType>& conditions,
                                                           StateType& kfwd,
                                                           StateType& dkfwd_dT,
                                                           VectorStateType& dkfwd_dX) const;

    // Deprecated API for backwards compatibility
    template <typename StateType, typename VectorStateType>
    void compute_forward_rate_coefficient_and_derivatives( const VectorStateType& molar_densities,
                                                           const
                                                           StateType& temp,
                                                           StateType& kfwd,
                                                           StateType& dkfwd_dT,
                                                           VectorStateType& dkfwd_dX) const;

    ////
    template <typename StateType, typename VectorStateType>
    StateType compute_rate_of_progress( const VectorStateType& molar_densities,
                                        const KineticsConditions<StateType,VectorStateType>& conditions,
                                        const StateType& P0_RT,
                                        const VectorStateType& h_RT_minus_s_R) const;

    // Deprecated API for backwards compatibility
    template <typename StateType, typename VectorStateType>
    StateType compute_rate_of_progress( const VectorStateType& molar_densities,
                                        const StateType& temp,
                                        const StateType& P0_RT,
                                        const VectorStateType& h_RT_minus_s_R) const;

    template <typename StateType, typename VectorStateType>
    void compute_rate_of_progress_and_derivatives( const VectorStateType &molar_densities,
                                                   const ChemicalMixture<CoeffType>& /*chem_mixture*/, // fully useless, why is it here?
                                                   const KineticsConditions<StateType,VectorStateType>& conditions,
                                                   const StateType &P0_RT,
                                                   const VectorStateType &h_RT_minus_s_R,
                                                   const VectorStateType &dh_RT_minus_s_R_dT,
                                                   StateType& net_reaction_rate,
                                                   StateType& dnet_rate_dT,
                                                   VectorStateType& dnet_rate_dX_s ) const;

    // Deprecated API for backwards compatibility
    template <typename StateType, typename VectorStateType>
    void compute_rate_of_progress_and_derivatives( const VectorStateType &molar_densities,
                                                   const ChemicalMixture<CoeffType>& chem_mixture,
                                                   const StateType& temp,
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

    //! Return writeable reference to the falloff object, test type
    //
    // just a wrapper around the FalloffType & F() method
    // of the FalloffReaction object, test consistency of type
    template <typename FalloffType>
    FalloffType & falloff();

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
    std::string _id;
    std::string _equation;
    std::vector<std::string> _reactant_names;
    std::vector<std::string> _product_names;
    std::vector<unsigned int> _reactant_ids;
    std::vector<unsigned int> _product_ids;
    std::vector<unsigned int> _reactant_stoichiometry;
    std::vector<unsigned int> _product_stoichiometry;
    std::vector<unsigned int> _species_reactant_stoichiometry;
    std::vector<unsigned int> _species_product_stoichiometry;
    std::vector<CoeffType>    _species_reactant_partial_order;
    std::vector<CoeffType>    _species_product_partial_order;
    std::vector<int>          _species_delta_stoichiometry;
    int _gamma;
    bool _initialized;
    bool _reversible;
    ReactionType::ReactionType _type;
    KineticsModel::KineticsModel _kintype;

    //! The forward reaction rate modified Arrhenius form.
    std::vector<KineticsType<CoeffType,VectorCoeffType>* > _forward_rate;

    //! efficiencies for three body reactions
    std::vector<CoeffType> _efficiencies;

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
  const std::string & Reaction<CoeffType,VectorCoeffType>::id() const
  {
    return _id;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  void Reaction<CoeffType,VectorCoeffType>::set_id(const std::string & id)
  {
    _id = id;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  const std::string & Reaction<CoeffType,VectorCoeffType>::equation() const
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
  CoeffType Reaction<CoeffType,VectorCoeffType>::reactant_partial_order(const unsigned int r) const
  {
    antioch_assert_less(r, _species_reactant_partial_order.size());
    antioch_assert_less(_reactant_ids[r], this->n_species());
    return _species_reactant_partial_order[r];
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  CoeffType Reaction<CoeffType,VectorCoeffType>::product_partial_order(const unsigned int p) const
  {
    antioch_assert_less(p, _species_product_partial_order.size());
    antioch_assert_less(_product_ids[p], this->n_species());
    return _species_product_partial_order[p];
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  void Reaction<CoeffType,VectorCoeffType>::add_reactant (const std::string &name,
                                          const unsigned int r_id,
                                          const unsigned int stoichiometric_coeff,
                                          const CoeffType partial_order)
  {
    antioch_assert_less(r_id, this->n_species());
    _reactant_names.push_back(name);
    _reactant_ids.push_back(r_id);
    _reactant_stoichiometry.push_back(stoichiometric_coeff);

   CoeffType order = (partial_order == std::numeric_limits<CoeffType>::infinity() )?static_cast<CoeffType>(stoichiometric_coeff):partial_order;
    _species_reactant_partial_order.push_back(order);
    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  void Reaction<CoeffType,VectorCoeffType>::add_product (const std::string &name,
                                         const unsigned int p_id,
                                         const unsigned int stoichiometric_coeff,
                                         const CoeffType partial_order)
  {
    antioch_assert_less(p_id, this->n_species());
    _product_names.push_back(name);
    _product_ids.push_back(p_id);
    _product_stoichiometry.push_back(stoichiometric_coeff);

   CoeffType order = (partial_order == std::numeric_limits<CoeffType>::infinity() )?static_cast<CoeffType>(stoichiometric_coeff):partial_order;
    _species_product_partial_order.push_back(order);
    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  void Reaction<CoeffType,VectorCoeffType>::clear_reactant()
  {
    _reactant_names.clear();
    _reactant_ids.clear();
    _reactant_stoichiometry.clear();
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  void Reaction<CoeffType,VectorCoeffType>::clear_product()
  {
    _product_names.clear();
    _product_ids.clear();
    _product_stoichiometry.clear();
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  void Reaction<CoeffType,VectorCoeffType>::set_efficiency (const std::string & , //? where does that come from?
                                            const unsigned int s,
                                            const CoeffType efficiency)
  {
    antioch_assert_less(s, _efficiencies.size());
    _efficiencies[s] = efficiency;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  CoeffType Reaction<CoeffType,VectorCoeffType>::get_efficiency (const unsigned int s) const
  {
    antioch_assert_less(s, _efficiencies.size());
    return _efficiencies[s];
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  CoeffType Reaction<CoeffType,VectorCoeffType>::efficiency( const unsigned int s ) const
  {

    antioch_assert(s < _efficiencies.size());
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
     return;
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
    _gamma = 0;
    for (unsigned int r=0; r< this->n_reactants(); r++)
      {
        _gamma -= this->reactant_stoichiometric_coefficient(r);
      }

    for (unsigned int p=0; p < this->n_products(); p++)
      {
        _gamma += this->product_stoichiometric_coefficient(p);
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
    antioch_assert( this->initialized() );
    //!\todo Make this assertion vector-compatible
    // antioch_assert_greater( P0_RT, 0.0 );
    antioch_assert_greater( h_RT_minus_s_R.size(), 0 );
    antioch_assert_equal_to( h_RT_minus_s_R.size(), this->n_species() );


// DrG0 = - reactants + product
// K = (P0/(RT))^gamma exp(-DrG0)
// exppower = -DrG0 = reactants - products
    StateType exppower = ( static_cast<CoeffType>(_reactant_stoichiometry[0])*
                            h_RT_minus_s_R[_reactant_ids[0]] );

    for (unsigned int s=1; s < this->n_reactants(); s++)
      {
        exppower += ( static_cast<CoeffType>(_reactant_stoichiometry[s])*
                       h_RT_minus_s_R[_reactant_ids[s]] );
      }

    for (unsigned int s=0; s < this->n_products(); s++)
      {
        exppower -= ( static_cast<CoeffType>(_product_stoichiometry[s])*
                       h_RT_minus_s_R[_product_ids[s]] );
      }
    return ant_pow( P0_RT, static_cast<CoeffType>(this->gamma()) ) * ant_exp(exppower);
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

    // get the equilibrium constant
    keq = this->equilibrium_constant( P0_RT, h_RT_minus_s_R );

    StateType ddT_exppower = ( static_cast<CoeffType>(_reactant_stoichiometry[0])*
                                ddT_h_RT_minus_s_R[_reactant_ids[0]] );

    for (unsigned int s=1; s<this->n_reactants(); s++)
      ddT_exppower +=  ( static_cast<CoeffType>(_reactant_stoichiometry[s])*
                         ddT_h_RT_minus_s_R[_reactant_ids[s]] );

    for (unsigned int s=0; s<this->n_products(); s++)
      ddT_exppower -=  ( static_cast<CoeffType>(_product_stoichiometry[s])*
                         ddT_h_RT_minus_s_R[_product_ids[s]] );

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
             << this->reactant_stoichiometric_coefficient(r) << ","
             << this->reactant_partial_order(r) << " ";
        os << "\n"
           << "#   products:  ";
        for (unsigned int p=0; p<this->n_products(); p++)
          os << this->product_name(p) << ":"
             << this->product_stoichiometric_coefficient(p) << ","
             << this->product_partial_order(p) << " ";
      }
    os << "\n#   Chemical process: " << _type;
    os << "\n#   Kinetics model: "   << _kintype;
    (_reversible)?os << "\n#   reversible":
                  os << "\n#   irreversible";
    for(unsigned int ir = 0; ir < _forward_rate.size(); ir++)
      {
        os << "\n#   forward rate eqn: " << *_forward_rate[ir];
      }

    if (_type == ReactionType::THREE_BODY ||
        _type == ReactionType::LINDEMANN_FALLOFF_THREE_BODY ||
        _type == ReactionType::TROE_FALLOFF_THREE_BODY)
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

      case(ReactionType::LINDEMANN_FALLOFF_THREE_BODY):
        {
          return (static_cast<const FalloffThreeBodyReaction<CoeffType,LindemannFalloff<CoeffType> >*>(this))->compute_forward_rate_coefficient(molar_densities,conditions);
        }
        break;

      case(ReactionType::TROE_FALLOFF_THREE_BODY):
        {
          return (static_cast<const FalloffThreeBodyReaction<CoeffType,TroeFalloff<CoeffType> >*>(this))->compute_forward_rate_coefficient(molar_densities,conditions);
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
  StateType Reaction<CoeffType,VectorCoeffType>::compute_forward_rate_coefficient
  ( const VectorStateType& molar_densities,
    const StateType& temp) const
  {
    antioch_deprecated();
    return compute_forward_rate_coefficient
      (molar_densities,
       KineticsConditions<StateType,VectorStateType>(temp));
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

      case(ReactionType::LINDEMANN_FALLOFF_THREE_BODY):
        {
          (static_cast<const FalloffThreeBodyReaction<CoeffType,LindemannFalloff<CoeffType> >*>(this))->compute_forward_rate_coefficient_and_derivatives(molar_densities,conditions,kfwd,dkfwd_dT,dkfwd_dX);
        }
        break;

      case(ReactionType::TROE_FALLOFF_THREE_BODY):
        {
          (static_cast<const FalloffThreeBodyReaction<CoeffType,TroeFalloff<CoeffType> >*>(this))->compute_forward_rate_coefficient_and_derivatives(molar_densities,conditions,kfwd,dkfwd_dT,dkfwd_dX);
        }
        break;

      default:
        {
          antioch_error();
        }

      } // switch(type)

    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  void
  Reaction<CoeffType,VectorCoeffType>::compute_forward_rate_coefficient_and_derivatives
  ( const VectorStateType& molar_densities,
    const StateType& temp,
    StateType& kfwd,
    StateType& dkfwd_dT,
    VectorStateType& dkfwd_dX) const
  {
    antioch_deprecated();
    return compute_forward_rate_coefficient_and_derivatives
      (molar_densities,
       KineticsConditions<StateType,VectorStateType>(temp),
       kfwd, dkfwd_dT, dkfwd_dX);
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
            this->reactant_partial_order(ro));
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
              this->product_partial_order(po));

        }
    }

    return kfwd_times_reactants - kbkwd_times_products;

  }

  template<typename CoeffType, typename VectorCoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  StateType Reaction<CoeffType,VectorCoeffType>::compute_rate_of_progress( const VectorStateType& molar_densities,
                                                           const StateType& temp,
                                                           const StateType& P0_RT,
                                                           const VectorStateType& h_RT_minus_s_R) const
  {
    antioch_deprecated();
    return compute_rate_of_progress
      (molar_densities,
       KineticsConditions<StateType,VectorStateType>(temp),
       P0_RT, h_RT_minus_s_R);
  }

  template<typename CoeffType, typename VectorCoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  void Reaction<CoeffType,VectorCoeffType>::compute_rate_of_progress_and_derivatives( const VectorStateType &molar_densities,
                                                                      const ChemicalMixture<CoeffType>& /*chem_mixture*/,
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
            this->reactant_partial_order(ro));

        const StateType dval =
          ( static_cast<CoeffType>(this->reactant_stoichiometric_coefficient(ro))*
            ant_pow( molar_densities[this->reactant_id(ro)],
              this->reactant_partial_order(ro) - 1)
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
              this->product_partial_order(po));

          const StateType dval =
            ( static_cast<CoeffType>(this->product_stoichiometric_coefficient(po))*
              ant_pow( molar_densities[this->product_id(po)],
              this->product_partial_order(po) - 1)
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

  template<typename CoeffType, typename VectorCoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  void Reaction<CoeffType,VectorCoeffType>::compute_rate_of_progress_and_derivatives( const VectorStateType &molar_densities,
                                                                      const ChemicalMixture<CoeffType>& chem_mixture,
                                                                      const StateType& temp,
                                                                      const StateType &P0_RT,
                                                                      const VectorStateType &h_RT_minus_s_R,
                                                                      const VectorStateType &dh_RT_minus_s_R_dT,
                                                                      StateType& net_reaction_rate,
                                                                      StateType& dnet_rate_dT,
                                                                      VectorStateType& dnet_rate_dX_s ) const
  {
    antioch_deprecated();
    return compute_rate_of_progress_and_derivatives
      (molar_densities, chem_mixture,
       KineticsConditions<StateType,VectorStateType>(temp),
       P0_RT, h_RT_minus_s_R, dh_RT_minus_s_R_dT, net_reaction_rate,
       dnet_rate_dT, dnet_rate_dX_s);
  }

  template<typename CoeffType, typename VectorCoeffType>
  template <typename FalloffType>
  inline
  FalloffType & Reaction<CoeffType,VectorCoeffType>::falloff()
  {
      switch(_type)
      {
        case(ReactionType::LINDEMANN_FALLOFF):
        {
          return (static_cast<const FalloffReaction<CoeffType,LindemannFalloff<CoeffType> >*>(this))->F();
        }
        break;

        case(ReactionType::TROE_FALLOFF):
        {
          return (static_cast<const FalloffReaction<CoeffType,LindemannFalloff<CoeffType> >*>(this))->F();
        }
        break;

        case(ReactionType::LINDEMANN_FALLOFF_THREE_BODY):
        {
          return (static_cast<const FalloffThreeBodyReaction<CoeffType,LindemannFalloff<CoeffType> >*>(this))->F();
        }
        break;

        case(ReactionType::TROE_FALLOFF_THREE_BODY):
        {
          return (static_cast<const FalloffThreeBodyReaction<CoeffType,LindemannFalloff<CoeffType> >*>(this))->F();
        }
        break;

        default:
        {
          std::cerr << "You are trying to retrieve a Falloff object in a reaction that is not a falloff.\n"
                    << "The reaction is " << *this << std::endl;
          antioch_error();
          return NULL;
        }

      }
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  void Reaction<CoeffType,VectorCoeffType>::set_parameter_of_rate(KineticsModel::Parameters parameter, CoeffType new_value, unsigned int n_kin, const std::string & unit)
  {
      antioch_assert_less(n_kin,_forward_rate.size());
      reset_parameter_of_rate(*_forward_rate[n_kin], parameter, new_value, unit);
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  void Reaction<CoeffType,VectorCoeffType>::set_parameter_of_rate(KineticsModel::Parameters parameter, CoeffType new_value, unsigned int n_kin, int l, const std::string & unit)
  {
      antioch_assert_less(n_kin,_forward_rate.size());
      reset_parameter_of_rate(*_forward_rate[n_kin], parameter, new_value, l, unit);
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  CoeffType Reaction<CoeffType,VectorCoeffType>::get_parameter_of_rate(KineticsModel::Parameters parameter, unsigned int n_reaction, const std::string & /* unit*/) const
  {
      antioch_assert_less(n_reaction,_forward_rate.size());
      return _forward_rate[n_reaction]->get_parameter(parameter);
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  CoeffType Reaction<CoeffType,VectorCoeffType>::get_parameter_of_rate(KineticsModel::Parameters parameter, unsigned int n_reaction, const std::string & /* unit*/, int l) const
  {
      antioch_assert_less(n_reaction,_forward_rate.size());
      return _forward_rate[n_reaction]->get_parameter(parameter,l);
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  void Reaction<CoeffType,VectorCoeffType>::set_parameter_of_chemical_process(ReactionType::Parameters parameter, CoeffType new_value, unsigned int species)
  {
      switch(parameter)
      {
        case ReactionType::EFFICIENCIES:
        {
          this->set_efficiency(std::string(),species,new_value); // tests inside, if not three-body or species wrong
        }
          break;
        case ReactionType::TROE_ALPHA:
        {
          if(this->type() == ReactionType::TROE_FALLOFF)
          {
             (static_cast<FalloffReaction<CoeffType,TroeFalloff<CoeffType> >*>(this))->F().set_alpha(new_value);
          }
          else if(this->type() == ReactionType::TROE_FALLOFF_THREE_BODY)
          {
             (static_cast<FalloffThreeBodyReaction<CoeffType,TroeFalloff<CoeffType> >*>(this))->F().set_alpha(new_value);
          }else
          {
            antioch_error();
          }
        }
          break;
        case ReactionType::TROE_T1:
        {
          if(this->type() == ReactionType::TROE_FALLOFF)
          {
             (static_cast<FalloffReaction<CoeffType,TroeFalloff<CoeffType> >*>(this))->F().set_T1(new_value);
          }
          else if(this->type() == ReactionType::TROE_FALLOFF_THREE_BODY)
          {
             (static_cast<FalloffThreeBodyReaction<CoeffType,TroeFalloff<CoeffType> >*>(this))->F().set_T1(new_value);
          }else
          {
            antioch_error();
          }
        }
          break;
        case ReactionType::TROE_T2:
        {
          if(this->type() == ReactionType::TROE_FALLOFF)
          {
             (static_cast<FalloffReaction<CoeffType,TroeFalloff<CoeffType> >*>(this))->F().set_T2(new_value);
          }
          else if(this->type() == ReactionType::TROE_FALLOFF_THREE_BODY)
          {
             (static_cast<FalloffThreeBodyReaction<CoeffType,TroeFalloff<CoeffType> >*>(this))->F().set_T2(new_value);
          }else
          {
            antioch_error();
          }
        }
          break;
        case ReactionType::TROE_T3:
        {
          if(this->type() == ReactionType::TROE_FALLOFF)
          {
             (static_cast<FalloffReaction<CoeffType,TroeFalloff<CoeffType> >*>(this))->F().set_T3(new_value);
          }
          else if(this->type() == ReactionType::TROE_FALLOFF_THREE_BODY)
          {
             (static_cast<FalloffThreeBodyReaction<CoeffType,TroeFalloff<CoeffType> >*>(this))->F().set_T3(new_value);
          }else
          {
            antioch_error();
          }
        }
          break;
        default:
        {
          antioch_error();
        }
          break;
      }
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  CoeffType Reaction<CoeffType,VectorCoeffType>::get_parameter_of_chemical_process(ReactionType::Parameters parameter, unsigned int species) const
  {
      CoeffType value(0);
      switch(parameter)
      {
        case ReactionType::EFFICIENCIES:
        {
          value = this->get_efficiency(species); // tests inside, if not three-body or species wrong
        }
          break;
        case ReactionType::TROE_ALPHA:
        {
          if(this->type() == ReactionType::TROE_FALLOFF)
          {
             value = (static_cast<const FalloffReaction<CoeffType,TroeFalloff<CoeffType> >*>(this))->F().get_alpha();
          }
          else if(this->type() == ReactionType::TROE_FALLOFF_THREE_BODY)
          {
             value = (static_cast<const FalloffThreeBodyReaction<CoeffType,TroeFalloff<CoeffType> >*>(this))->F().get_alpha();
          }else
          {
            antioch_error();
          }
        }
          break;
        case ReactionType::TROE_T1:
        {
          if(this->type() == ReactionType::TROE_FALLOFF)
          {
             value = (static_cast<const FalloffReaction<CoeffType,TroeFalloff<CoeffType> >*>(this))->F().get_T1();
          }
          else if(this->type() == ReactionType::TROE_FALLOFF_THREE_BODY)
          {
             value = (static_cast<const FalloffThreeBodyReaction<CoeffType,TroeFalloff<CoeffType> >*>(this))->F().get_T1();
          }else
          {
            antioch_error();
          }
        }
          break;
        case ReactionType::TROE_T2:
        {
          if(this->type() == ReactionType::TROE_FALLOFF)
          {
             value = (static_cast<const FalloffReaction<CoeffType,TroeFalloff<CoeffType> >*>(this))->F().get_T2();
          }
          else if(this->type() == ReactionType::TROE_FALLOFF_THREE_BODY)
          {
             value = (static_cast<const FalloffThreeBodyReaction<CoeffType,TroeFalloff<CoeffType> >*>(this))->F().get_T2();
          }else
          {
            antioch_error();
          }
        }
          break;
        case ReactionType::TROE_T3:
        {
          if(this->type() == ReactionType::TROE_FALLOFF)
          {
             value = (static_cast<const FalloffReaction<CoeffType,TroeFalloff<CoeffType> >*>(this))->F().get_T3();
          }
          else if(this->type() == ReactionType::TROE_FALLOFF_THREE_BODY)
          {
             value = (static_cast<const FalloffThreeBodyReaction<CoeffType,TroeFalloff<CoeffType> >*>(this))->F().get_T3();
          }else
          {
            antioch_error();
          }
        }
          break;
        default:
        {
          antioch_error();
        }
          break;
      }

    return value;
  }

} // namespace Antioch

#endif // ANTIOCH_REACTION_H
