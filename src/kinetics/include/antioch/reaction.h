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

//C++
#include <string>
#include <vector>
#include <iostream>

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/reaction_enum.h"
#include "antioch/arrhenius_rate.h"

namespace Antioch
{
  //!A single reaction mechanism. 
  /*!
    This class encapsulates a single reaction mechanism.  The mechanism could be 
    an elementary reaction, or a three-body reaction.  All reactions are assumed
    to be reversible. This class was originally taken from \p FIN-S.
    \todo{Do we want to template this class around the rate type?}
  */
  template<class CoeffType>
  class Reaction
  {
  public:

    //! Construct a single reaction mechanism.
    Reaction( const unsigned int n_species, const std::string &equation );
    
    ~Reaction();

    unsigned int n_species() const;
    
    //! \returns the equation for this reaction.
    std::string equation() const;
    
    //! Type of reaction. 
    /*! Type of reaction. Presently only ELEMENTARY or THREE_BODY
     *  reversible reactions are considered.
     */
    ReactionType::ReactionType type() const;

    //! Set the type of reaction. 
    /*! Set the type of reaction. Presently only ELEMENTARY or THREE_BODY
     * reversible reactions are considered.
     */
    void set_type( const ReactionType::ReactionType type);

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
    template <typename StateType>
    StateType equilibrium_constant( const StateType P0_RT,
				    const std::vector<StateType>& h_RT_minus_s_R ) const;

    //!
    template <typename StateType>
    void equilibrium_constant_and_derivative( const StateType T,
					      const StateType P0_RT,
					      const std::vector<StateType>& h_RT_minus_s_R,
					      const std::vector<StateType>& ddT_h_RT_minus_s_R,
					      StateType& keq,
					      StateType& dkeq_dT) const;

    //!
    template <typename StateType>
    StateType compute_rate_of_progress( const std::vector<StateType>& molar_densities,
					const StateType kfwd, 
					const StateType kbkwd ) const;
    
    //!
    template <typename StateType>
    void compute_rate_of_progress_and_derivatives( const std::vector<StateType>& molar_densities,
						   const std::vector<StateType>& molar_mass,
						   const StateType kfwd,  
						   const StateType dkfwd_dT, 
						   const StateType kbkwd,
						   const StateType dkbkwd_dT,
						   StateType& Rfwd,
						   StateType& dRfwd_dT,
						   std::vector<StateType>& dRfwd_drho, 
						   StateType& Rbkwd,
						   StateType& dRbkwd_dT,
						   std::vector<StateType>& dRbkwd_drho) const;

    //! Return const reference to the forward rate object
    const ArrheniusRate<CoeffType>& forward_rate() const;

    //! Return writeable reference to the forward rate object
    ArrheniusRate<CoeffType>& forward_rate();

    //! Formatted print, by default to \p std::cout.
    void print(std::ostream& os = std::cout) const;

    //! Formatted print.
    friend std::ostream& operator << (std::ostream& os, const Reaction &rxn)
    {
      rxn.print(os);
      return os;
    }

  private:
    
    unsigned int _n_species;
    ReactionType::ReactionType _type;
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

    //! The forward reaction rate modified Arrhenius form.
    ArrheniusRate<CoeffType> _forward_rate;

  };

  /* ------------------------- Inline Functions -------------------------*/
  template<class CoeffType>
  inline
  unsigned int Reaction<CoeffType>::n_species() const
  {
    return _n_species;
  }

  template<class CoeffType>
  inline
  std::string Reaction<CoeffType>::equation() const
  {
    return _equation;
  }

  template<class CoeffType>
  inline
  ReactionType::ReactionType Reaction<CoeffType>::type() const
  {
    return _type;
  }

  template<class CoeffType>
  inline
  void Reaction<CoeffType>::set_type( const ReactionType::ReactionType type)
  {
    _type = type;
    return;
  }

  template<class CoeffType>
  inline
  bool Reaction<CoeffType>::initialized() const
  {
    return _initialized;
  }

  template<class CoeffType>
  inline
  unsigned int Reaction<CoeffType>::n_reactants () const
  {
    antioch_assert_less(_reactant_ids.size(), this->n_species());
    antioch_assert_equal_to(_reactant_ids.size(), _reactant_stoichiometry.size());
    antioch_assert_equal_to(_reactant_ids.size(), _reactant_names.size());
    return _reactant_ids.size();
  }

  template<class CoeffType>
  inline
  unsigned int Reaction<CoeffType>::n_products() const
  {
    antioch_assert_less(_product_ids.size(), this->n_species());
    antioch_assert_equal_to(_product_ids.size(), _product_stoichiometry.size());
    antioch_assert_equal_to(_product_ids.size(),  _product_names.size());
    return _product_ids.size();
  }

  template<class CoeffType>
  inline
  const std::string& Reaction<CoeffType>::reactant_name(const unsigned int r) const
  {       
    antioch_assert_less(r, _reactant_names.size()); 
    return _reactant_names[r]; 
  }

  template<class CoeffType>
  inline
  const std::string& Reaction<CoeffType>::product_name(const unsigned int p) const
  { 
    antioch_assert_less(p, _product_names.size()); 
    return _product_names[p]; 
  }
  
  template<class CoeffType>
  inline
  unsigned int Reaction<CoeffType>::reactant_id(const unsigned int r) const
  { 
    antioch_assert_less(r, _reactant_ids.size()); 
    antioch_assert_less(_reactant_ids[r], this->n_species());
    return _reactant_ids[r]; 
  }

  template<class CoeffType>
  inline
  unsigned int Reaction<CoeffType>::product_id(const unsigned int p) const
  { 
    antioch_assert_less(p, _product_ids.size()); 
    antioch_assert_less(_product_ids[p], this->n_species());
    return _product_ids[p]; 
  }

  template<class CoeffType>
  inline
  unsigned int Reaction<CoeffType>::reactant_stoichiometric_coefficient(const unsigned int r) const
  {      
    antioch_assert_less(r, _reactant_stoichiometry.size());
    antioch_assert_less(_reactant_ids[r], this->n_species());
    return _reactant_stoichiometry[r];
  }

  template<class CoeffType>
  inline
  unsigned int Reaction<CoeffType>::product_stoichiometric_coefficient(const unsigned int p) const
  {
    antioch_assert_less(p, _product_stoichiometry.size());
    antioch_assert_less(_product_ids[p], this->n_species());
    return _product_stoichiometry[p];
  }

  template<class CoeffType>
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

  template<class CoeffType>
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

  template<class CoeffType>
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

  template<class CoeffType>
  inline
  CoeffType Reaction<CoeffType>::efficiency( const unsigned int s ) const
  {
    antioch_assert_less(s, _efficiencies.size());
    antioch_assert_equal_to(_type, ReactionType::THREE_BODY);
    return _efficiencies[s];
  }

  template<class CoeffType>
  inline
  int Reaction<CoeffType>::gamma () const
  {
    return _gamma;
  }

  template<class CoeffType>
  inline
  const ArrheniusRate<CoeffType>& Reaction<CoeffType>::forward_rate() const
  {
    return _forward_rate;
  }

  template<class CoeffType>
  inline
  ArrheniusRate<CoeffType>& Reaction<CoeffType>::forward_rate()
  {
    return _forward_rate;
  }


  template<class CoeffType>
  inline
  Reaction<CoeffType>::Reaction( const unsigned int n_species,
				 const std::string &equation ) 
    : _n_species(n_species),
      _type(ReactionType::ELEMENTARY),
      _equation(equation),
      _gamma(0),
      _initialized(false)
  {
    _efficiencies.resize(_n_species); 
    std::fill (_efficiencies.begin(), _efficiencies.end(), 1.);
  }


  template<class CoeffType>
  inline
  Reaction<CoeffType>::~Reaction()
  {
    return;
  }


  template<class CoeffType>
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


  template<class CoeffType>
  template<class StateType>
  inline
  StateType Reaction<CoeffType>::equilibrium_constant( const StateType P0_RT,
						       const std::vector<StateType>& h_RT_minus_s_R ) const
  {
    antioch_assert( this->initialized() );
    antioch_assert_greater( P0_RT, 0.0 );
    antioch_assert_equal_to( h_RT_minus_s_R.size(), this->n_species() );
    antioch_assert_equal_to( _species_delta_stoichiometry.size(), this->n_species() );

    StateType exppower = 0.;

    for (unsigned int s=0; s < this->n_species(); s++)
      {
	exppower += -( static_cast<CoeffType>(_species_delta_stoichiometry[s])*
		       h_RT_minus_s_R[s] );
      }

    return std::pow( P0_RT, this->gamma() )*std::exp(exppower);
  }


  template<class CoeffType>
  template<class StateType>
  inline
  void Reaction<CoeffType>::equilibrium_constant_and_derivative( const StateType T,
								 const StateType P0_RT,
								 const std::vector<StateType>& h_RT_minus_s_R,
								 const std::vector<StateType>& ddT_h_RT_minus_s_R,
								 StateType& keq,
								 StateType& dkeq_dT) const
  {
    antioch_assert(this->initialized());
    antioch_assert_greater( P0_RT, 0.0 );
    antioch_assert_greater( T, 0.0 );
    antioch_assert_equal_to( h_RT_minus_s_R.size(), this->n_species() );
    antioch_assert_equal_to( ddT_h_RT_minus_s_R.size(), this->n_species() );
    antioch_assert_equal_to( _species_delta_stoichiometry.size(), this->n_species() );

    // get the equilibrium constant
    keq = this->equilibrium_constant( P0_RT, h_RT_minus_s_R );

    StateType ddT_exppower = 0.;

    for (unsigned int s=0; s<this->n_species(); s++)
      ddT_exppower += -( static_cast<CoeffType>(_species_delta_stoichiometry[s])*
			 ddT_h_RT_minus_s_R[s] );

    // compute its derivative
    dkeq_dT = keq*(-static_cast<CoeffType>(this->gamma())/T + ddT_exppower);

    return;
  }


  template<class CoeffType>
  template<class StateType>
  inline
  StateType Reaction<CoeffType>::compute_rate_of_progress( const std::vector<StateType>& molar_densities,
							   const StateType kfwd, 
							   const StateType kbkwd ) const
  {
    antioch_assert_equal_to( molar_densities.size(), this->n_species() );

    StateType Rfwd = 0.0;
    StateType Rbkwd = 0.0;
    StateType kfwd_times_reactants = kfwd;
    StateType kbkwd_times_products = kbkwd;

    for (unsigned int r=0; r<this->n_reactants(); r++)
      {
	kfwd_times_reactants *= std::pow( molar_densities[this->reactant_id(r)],
					  static_cast<int>(this->reactant_stoichiometric_coefficient(r)) );
      }
      
    for (unsigned int p=0; p<this->n_products(); p++)
      {
	kbkwd_times_products *= std::pow( molar_densities[this->product_id(p)],
					  static_cast<int>(this->product_stoichiometric_coefficient(p)) );
      }
      
    switch (this->type())
      {
	// for elementary reactions the forward and backward
	// rates of progress are simply the rates times
	// the product of the reactants, products respectively.
      case(ReactionType::ELEMENTARY):
	{
	  Rfwd  = kfwd_times_reactants;
	  Rbkwd = kbkwd_times_products;
	}
	break;
	
	// for threebody reactions we need to include the
	// contrbution from each collision partner
      case(ReactionType::THREE_BODY):
	{
	  for (unsigned int s=0; s<this->n_species(); s++)
	    {	     
	      Rfwd += ( kfwd_times_reactants * 
			this->efficiency(s) *
			molar_densities[s] );

	      Rbkwd += ( kbkwd_times_products *
			 this->efficiency(s) *
			 molar_densities[s] );
	    }
	}
	break;
	
      default:
	{
	  std::cerr << "Error: Invalid reaction type " << this->type() << std::endl;
	  antioch_error();
	}
	break;
      }

    // Rnet = Rfwd - Rbkwd
    return (Rfwd - Rbkwd);
  }


  template<class CoeffType>
  template<class StateType>
  inline
  void Reaction<CoeffType>::compute_rate_of_progress_and_derivatives( const std::vector<StateType> &molar_densities,
								      const std::vector<StateType> &molar_mass,
								      const StateType kfwd, 
								      const StateType dkfwd_dT,
								      const StateType kbkwd,
								      const StateType dkbkwd_dT,
								      StateType& Rfwd,
								      StateType& dRfwd_dT,
								      std::vector<StateType>& dRfwd_drho, 
								      StateType& Rbkwd,
								      StateType& dRbkwd_dT, std::vector<StateType> &dRbkwd_drho) const
  {
    antioch_assert_equal_to (molar_densities.size(), this->n_species());
    antioch_assert_equal_to (molar_mass.size(),      this->n_species());

    Rfwd = 0.0;
    dRfwd_dT = 0.0;
    Rbkwd = 0.0;
    dRbkwd_dT = 0.0;
    dRfwd_drho.resize(this->n_species());
    dRbkwd_drho.resize (this->n_species());

    std::fill( dRfwd_drho.begin(),  dRfwd_drho.end(),  0.);    
    std::fill( dRbkwd_drho.begin(), dRbkwd_drho.end(), 0.);
    
    StateType kfwd_times_reactants = kfwd;
    StateType kbkwd_times_products = kbkwd;
    StateType ddT_kfwd_times_reactants = dkfwd_dT;
    StateType ddT_kbkwd_times_products = dkbkwd_dT;
      
    // pre-fill the participating species partials with the rates
    for (unsigned int r=0; r< this->n_reactants(); r++)
      {
	dRfwd_drho[this->reactant_id(r)] = kfwd;
      }
    
    for (unsigned int p=0; p < this->n_products(); p++)
      {
	dRbkwd_drho[this->product_id(p)] = kbkwd;
      }
    
    // Rfwd & derivatives
    for (unsigned int ro=0; ro < this->n_reactants(); ro++)
      {
	const StateType val = 
	  std::pow( molar_densities[this->reactant_id(ro)],
		    static_cast<int>(this->reactant_stoichiometric_coefficient(ro)) );
	  
	const StateType dval = 
	  ( static_cast<CoeffType>(this->reactant_stoichiometric_coefficient(ro))*
	    std::pow( molar_densities[this->reactant_id(ro)],
		      static_cast<int>(this->reactant_stoichiometric_coefficient(ro))-1 ) 
	    / molar_mass[this->reactant_id(ro)] );
	  	  
	kfwd_times_reactants     *= val;
	ddT_kfwd_times_reactants *= val;

	for (unsigned int ri=0; ri<this->n_reactants(); ri++)
	  {
	    dRfwd_drho[this->reactant_id(ri)] *= (ri == ro) ? dval : val;
	  }
      }

    // Rbkwd & derivatives
    for (unsigned int po=0; po< this->n_products(); po++)
      {
	const StateType val = 
	  std::pow( molar_densities[this->product_id(po)],
		    static_cast<int>(this->product_stoichiometric_coefficient(po)) );
	  
	const StateType dval = 
	  ( static_cast<CoeffType>(this->product_stoichiometric_coefficient(po))*
	    std::pow( molar_densities[this->product_id(po)],
		      static_cast<int>(this->product_stoichiometric_coefficient(po))-1 )
	    / molar_mass[this->product_id(po)] );
	
	kbkwd_times_products     *= val;
	ddT_kbkwd_times_products *= val;
	
	for (unsigned int pi=0; pi<this->n_products(); pi++)
	  {
	    dRbkwd_drho[this->product_id(pi)] *= (pi == po) ? dval : val;
	  }
      }

    switch (this->type())
      {
	// for elementary reactions the forward and backward
	// rates of progress are simply the rates times
	// the product of the reactants, products respectively.
      case(ReactionType::ELEMENTARY):
	{
	  Rfwd  = kfwd_times_reactants;
	  Rbkwd = kbkwd_times_products;

	  dRfwd_dT  = ddT_kfwd_times_reactants;
	  dRbkwd_dT = ddT_kbkwd_times_products;

	  // and the derivatives are already handled.
	}
	break;
	
	// for threebody reactions we need to include the
	// contrbution from each collision partner
      case(ReactionType::THREE_BODY):
	{
	  StateType summed_value=0.0;

	  for (unsigned int s=0; s < this->n_species(); s++)
	    {
	      summed_value += this->efficiency(s) * molar_densities[s];
	    }
	  
	  Rfwd = kfwd_times_reactants * summed_value;
	  
	  dRfwd_dT = ddT_kfwd_times_reactants * summed_value;
	  
	  Rbkwd = kbkwd_times_products * summed_value;
	  
	  dRbkwd_dT = ddT_kbkwd_times_products * summed_value;
	  
	  for (unsigned int s=0; s<this->n_species(); s++)
	    {
	      dRfwd_drho[s]  *= summed_value;
	      dRbkwd_drho[s] *= summed_value;
	      
	      // and the efficiency contribution derivative
	      dRfwd_drho[s]  += 
		this->efficiency(s) / molar_mass[s] * kfwd_times_reactants;

	      dRbkwd_drho[s] +=
		this->efficiency(s) / molar_mass[s] * kbkwd_times_products;
	    }
	}
	break;
	
      default:
	{
	  std::cerr << "Error: Invalid reaction type " << this->type() << std::endl;
	  antioch_error();
	}
	break;
      }

    return;
  }


  template<class CoeffType>
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
    os << "\n#   forward rate eqn: " << _forward_rate;
    
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
