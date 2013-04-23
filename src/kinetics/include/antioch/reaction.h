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

#ifndef ANTIOCH_REACTION_H
#define ANTIOCH_REACTION_H

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/kinetics_type.h"
#include "antioch/hercourtessen_rate.h"
#include "antioch/berthelot_rate.h"
#include "antioch/arrhenius_rate.h"
#include "antioch/berthelothercourtessen_rate.h"
#include "antioch/kooij_rate.h"
#include "antioch/vanthoff_rate.h"
#include "antioch/reaction_enum.h"

//C++
#include <string>
#include <vector>
#include <iostream>

namespace Antioch
{
  //!A single reaction mechanism. 
  /*!\class Reaction
 *
    This virtual base is derived for the following processes:
        - elementary process,
        - duplicate process,
        - falloff processes with:
                - Lindemann falloff,
                - Troe falloff.
    This class encapsulates a kinetics model.  The choosable kinetics models are
        - Hercourt Hessen \f$\alpha(T) = A T^\beta\f$
        - Berthelot \f$\alpha(T) = A \exp\left(D T\right)\f$
        - Arrhenius \f$\alpha(T) = A \exp\left(-\frac{E_a}{T}\right)\f$
        - Berthelot Hercourt Hessen \f$\alpha(T) = A T^\beta \exp\left(D T\right)\f$
        - Kooij \f$\alpha(T) = A T^\beta \exp\left(- \frac{E_a}{T}\right)\f$
        - Van't Hoff \f$\alpha(T) = A T^\beta \exp\left(- \frac{E_a}{T} + D T\right)\f$
    All reactions are assumed to be reversible. 
    By default, we choose an elementary process with a Kooij equation.
  */
  template<typename CoeffType=double>
  class Reaction
  {
  public:

    //! Construct a single reaction mechanism.
    Reaction( const unsigned int n_species, const std::string &equation, const ReactionType::ReactionType type = ReactionType::ELEMENTARY, const KinMod::KinMod kin = KinMod::KOOIJ);
    
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

    //! Model of kinetics.
    KinMod::KinMod kinetics_model() const;

    //! Set the model of kinetics.
    void set_kinetics_model( const KinMod::KinMod kin);
    
    bool initialized() const;

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
    void initialize();    

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


//// in reaction set

    //!
    template <typename StateType, typename VectorStateType>
    StateType compute_forward_rate_coefficient( const VectorStateType& molar_densities,
                                        const StateType& T) const;
    
    //!
    template <typename StateType, typename VectorStateType>
    void compute_forward_rate_coefficient_and_derivatives( const VectorStateType& molar_densities,
                                                           const StateType& T, 
                                                           StateType& kfwd, 
                                                           StateType& dkfwd_dT,
                                                           VectorStateType& dkfwd_dY) const;
////
    //! Return const reference to the forward rate object
    const KineticsType<CoeffType>& forward_rate(unsigned int ir = 0) const;

    //! Return writeable reference to the forward rate object
    KineticsType<CoeffType>& forward_rate(unsigned int ir = 0);

    //! Return const reference to the forward rate object
    void add_forward_rate(KineticsType<CoeffType> *rate);

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
    ReactionType::ReactionType _type;
    KinMod::KinMod _kintype;

    //! The forward reaction rate modified Arrhenius form.
    std::vector<KineticsType<CoeffType>* > _forward_rate;

  };

  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType>
  inline
  unsigned int Reaction<CoeffType>::n_species() const
  {
    return _n_species;
  }

  template<typename CoeffType>
  inline
  std::string Reaction<CoeffType>::equation() const
  {
    return _equation;
  }

  template<typename CoeffType>
  inline
  ReactionType::ReactionType Reaction<CoeffType>::type() const
  {
    return _type;
  }

  template<typename CoeffType>
  inline
  void Reaction<CoeffType>::set_type( const ReactionType::ReactionType type)
  {
    _type = type;
    return;
  }

  template<typename CoeffType>
  inline
  KinMod::KinMod Reaction<CoeffType>::kinetics_model() const
  {
    return _kintype;
  }

  template<typename CoeffType>
  inline
  void Reaction<CoeffType>::set_kinetics_model( const KinMod::KinMod kin)
  {
     _kintype = kin;
     return;
  }

  template<typename CoeffType>
  inline
  bool Reaction<CoeffType>::initialized() const
  {
    return _initialized;
  }


  template<typename CoeffType>
  inline
  void Reaction<CoeffType>::add_forward_rate(KineticsType<CoeffType> *rate)
  {
    _forward_rate.push_back(rate);
    return;
  }

  template<typename CoeffType>
  inline
  unsigned int Reaction<CoeffType>::n_reactants () const
  {
    antioch_assert_less(_reactant_ids.size(), this->n_species());
    antioch_assert_equal_to(_reactant_ids.size(), _reactant_stoichiometry.size());
    antioch_assert_equal_to(_reactant_ids.size(), _reactant_names.size());
    return _reactant_ids.size();
  }

  template<typename CoeffType>
  inline
  unsigned int Reaction<CoeffType>::n_products() const
  {
    antioch_assert_less(_product_ids.size(), this->n_species());
    antioch_assert_equal_to(_product_ids.size(), _product_stoichiometry.size());
    antioch_assert_equal_to(_product_ids.size(),  _product_names.size());
    return _product_ids.size();
  }

  template<typename CoeffType>
  inline
  const std::string& Reaction<CoeffType>::reactant_name(const unsigned int r) const
  {       
    antioch_assert_less(r, _reactant_names.size()); 
    return _reactant_names[r]; 
  }

  template<typename CoeffType>
  inline
  const std::string& Reaction<CoeffType>::product_name(const unsigned int p) const
  { 
    antioch_assert_less(p, _product_names.size()); 
    return _product_names[p]; 
  }
  
  template<typename CoeffType>
  inline
  unsigned int Reaction<CoeffType>::reactant_id(const unsigned int r) const
  { 
    antioch_assert_less(r, _reactant_ids.size()); 
    antioch_assert_less(_reactant_ids[r], this->n_species());
    return _reactant_ids[r]; 
  }

  template<typename CoeffType>
  inline
  unsigned int Reaction<CoeffType>::product_id(const unsigned int p) const
  { 
    antioch_assert_less(p, _product_ids.size()); 
    antioch_assert_less(_product_ids[p], this->n_species());
    return _product_ids[p]; 
  }

  template<typename CoeffType>
  inline
  unsigned int Reaction<CoeffType>::reactant_stoichiometric_coefficient(const unsigned int r) const
  {      
    antioch_assert_less(r, _reactant_stoichiometry.size());
    antioch_assert_less(_reactant_ids[r], this->n_species());
    return _reactant_stoichiometry[r];
  }

  template<typename CoeffType>
  inline
  unsigned int Reaction<CoeffType>::product_stoichiometric_coefficient(const unsigned int p) const
  {
    antioch_assert_less(p, _product_stoichiometry.size());
    antioch_assert_less(_product_ids[p], this->n_species());
    return _product_stoichiometry[p];
  }

  template<typename CoeffType>
  inline
  void Reaction<CoeffType>::add_reactant (const std::string &name,
                                          const unsigned int r_id,
                                          const unsigned int stoichiometric_coeff)
  {
    antioch_assert_less(r_id, this->n_species());
    _reactant_names.push_back(name);
    _reactant_ids.push_back(r_id);
    _reactant_stoichiometry.push_back(stoichiometric_coeff);
    return;
  }

  template<typename CoeffType>
  inline
  void Reaction<CoeffType>::add_product (const std::string &name,
                                         const unsigned int p_id,
                                         const unsigned int stoichiometric_coeff)
  {
    antioch_assert_less(p_id, this->n_species());
    _product_names.push_back(name);
    _product_ids.push_back(p_id);
    _product_stoichiometry.push_back(stoichiometric_coeff);
    return;
  }

  template<typename CoeffType>
  inline
  void Reaction<CoeffType>::set_efficiency (const std::string &,
                                            const unsigned int s,
                                            const CoeffType efficiency)
  {
    antioch_assert_less(s, this->n_species());
    antioch_assert_less(s, _efficiencies.size());
    antioch_assert_equal_to(_type, ReactionType::THREE_BODY);
    _efficiencies[s] = efficiency;
    return;
  }

  template<typename CoeffType>
  inline
  CoeffType Reaction<CoeffType>::efficiency( const unsigned int s ) const
  {
    antioch_assert_less(s, _efficiencies.size());
    antioch_assert_equal_to(_type, ReactionType::THREE_BODY);
    return _efficiencies[s];
  }

  template<typename CoeffType>
  inline
  int Reaction<CoeffType>::gamma () const
  {
    return _gamma;
  }

  template<typename CoeffType>
  inline
  const KineticsType<CoeffType>& Reaction<CoeffType>::forward_rate(unsigned int ir) const
  {
    return *_forward_rate[ir];
  }

  template<typename CoeffType>
  inline
  KineticsType<CoeffType>& Reaction<CoeffType>::forward_rate(unsigned int ir)
  {
    return *_forward_rate[ir];
  }


  template<typename CoeffType>
  inline
  Reaction<CoeffType>::Reaction( const unsigned int n_species,
                                 const std::string &equation, 
                                 const ReactionType::ReactionType type,
                                 const KinMod::KinMod kin) 
    : _n_species(n_species),
      _equation(equation),
      _gamma(0),
      _initialized(false),
      _type(type),
      _kintype(kin)
  {
    _efficiencies.resize(_n_species); 
    std::fill (_efficiencies.begin(), _efficiencies.end(), 1.);
  }


  template<typename CoeffType>
  inline
  Reaction<CoeffType>::~Reaction()
  {
    for(unsigned int ir = 0; ir < _forward_rate.size(); ir++)delete _forward_rate[ir];
    return;
  }


  template<typename CoeffType>
  inline
  void Reaction<CoeffType>::initialize()
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
    for (unsigned int s=0, _gamma=0; s<this->n_species(); s++)
      {
        _species_delta_stoichiometry[s] =
          ( _species_product_stoichiometry[s] -
            _species_reactant_stoichiometry[s] );
        
        _gamma += _species_delta_stoichiometry[s];
      }
   
    // set initialization flag
    _initialized = true;
  }


  template<typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  StateType Reaction<CoeffType>::equilibrium_constant( const StateType& P0_RT,
                                                       const VectorStateType& h_RT_minus_s_R ) const
  {
    using std::exp;
    using std::pow;

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

    return pow( P0_RT, static_cast<CoeffType>(this->gamma()) )*exp(exppower);
  }


  template<typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void Reaction<CoeffType>::equilibrium_constant_and_derivative( const StateType& T,
                                                                 const StateType& P0_RT,
                                                                 const VectorStateType& h_RT_minus_s_R,
                                                                 const VectorStateType& ddT_h_RT_minus_s_R,
                                                                 StateType& keq,
                                                                 StateType& dkeq_dT) const
  {
    antioch_assert(this->initialized());
    antioch_assert_greater( P0_RT, 0.0 );
    antioch_assert_greater( T, 0.0 );
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


  template<typename CoeffType>
  inline
  void Reaction<CoeffType>::print( std::ostream& os ) const
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
  
} // namespace Antioch

#endif // ANTIOCH_REACTION_H
