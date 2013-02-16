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

//This class
#include "antioch/reaction.h"

namespace Antioch
{
  template<class NumericType>
  Reaction<NumericType>::Reaction( const unsigned int n_species,
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

  template<class NumericType>
  Reaction<NumericType>::~Reaction()
  {
    return;
  }

  template<class NumericType>
  void Reaction<NumericType>::initialize()
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

  template<class NumericType>
  NumericType Reaction<NumericType>::equilibrium_constant( const NumericType P0_RT,
							   const std::vector<NumericType>& h_RT_minus_s_R ) const
  {
    antioch_assert( this->initialized() );
    antioch_assert_greater( P0_RT, 0.0 );
    antioch_assert_equal_to( h_RT_minus_s_R.size(), this->n_species() );
    antioch_assert_equal_to( _species_delta_stoichiometry.size(), this->n_species() );

    NumericType exppower = 0.;

    for (unsigned int s=0; s < this->n_species(); s++)
      {
	exppower += -( static_cast<NumericType>(_species_delta_stoichiometry[s])*
		       h_RT_minus_s_R[s] );
      }

    return std::pow( P0_RT, this->gamma() )*std::exp(exppower);
  }

  template<class NumericType>
  void Reaction<NumericType>::equilibrium_constant_and_derivative( const NumericType T,
								   const NumericType P0_RT,
								   const std::vector<NumericType>& h_RT_minus_s_R,
								   const std::vector<NumericType>& ddT_h_RT_minus_s_R,
								   NumericType& keq,
								   NumericType& dkeq_dT) const
  {
    antioch_assert(this->initialized());
    antioch_assert_greater( P0_RT, 0.0 );
    antioch_assert_greater( T, 0.0 );
    antioch_assert_equal_to( h_RT_minus_s_R.size(), this->n_species() );
    antioch_assert_equal_to( ddT_h_RT_minus_s_R.size(), this->n_species() );
    antioch_assert_equal_to( _species_delta_stoichiometry.size(), this->n_species() );

    // get the equilibrium constant
    keq = this->equilibrium_constant( P0_RT, h_RT_minus_s_R );

    NumericType ddT_exppower = 0.;

    for (unsigned int s=0; s<this->n_species(); s++)
      ddT_exppower += -( static_cast<NumericType>(_species_delta_stoichiometry[s])*
			 ddT_h_RT_minus_s_R[s] );

    // compute its derivative
    dkeq_dT = keq*(-static_cast<NumericType>(this->gamma())/T + ddT_exppower);

    return;
  }

  template<class NumericType>
  NumericType Reaction<NumericType>::compute_rate_of_progress( const std::vector<NumericType>& molar_densities,
							       const NumericType kfwd, 
							       const NumericType kbkwd ) const
  {
    antioch_assert_equal_to( molar_densities.size(), this->n_species() );

    NumericType Rfwd = 0.0;
    NumericType Rbkwd = 0.0;
    NumericType kfwd_times_reactants = kfwd;
    NumericType kbkwd_times_products = kbkwd;

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

  template<class NumericType>
  void Reaction<NumericType>::compute_rate_of_progress_and_derivatives( const std::vector<NumericType> &molar_densities,
									const std::vector<NumericType> &molar_mass,
									const NumericType kfwd, 
									const NumericType dkfwd_dT,
									const NumericType kbkwd,
									const NumericType dkbkwd_dT,
									NumericType& Rfwd,
									NumericType& dRfwd_dT,
									std::vector<NumericType>& dRfwd_drho, 
									NumericType& Rbkwd,
									NumericType& dRbkwd_dT, std::vector<NumericType> &dRbkwd_drho) const
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
    
    NumericType kfwd_times_reactants = kfwd;
    NumericType kbkwd_times_products = kbkwd;
    NumericType ddT_kfwd_times_reactants = dkfwd_dT;
    NumericType ddT_kbkwd_times_products = dkbkwd_dT;
      
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
	const NumericType val = 
	  std::pow( molar_densities[this->reactant_id(ro)],
		    static_cast<int>(this->reactant_stoichiometric_coefficient(ro)) );
	  
	const NumericType dval = 
	  ( static_cast<NumericType>(this->reactant_stoichiometric_coefficient(ro))*
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
	const NumericType val = 
	  std::pow( molar_densities[this->product_id(po)],
		    static_cast<int>(this->product_stoichiometric_coefficient(po)) );
	  
	const NumericType dval = 
	  ( static_cast<NumericType>(this->product_stoichiometric_coefficient(po))*
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
	  NumericType summed_value=0.0;

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

  template<class NumericType>
  void Reaction<NumericType>::print( std::ostream& os ) const
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
  
  /* ------------------------- Instantiate ------------------------- */
  template class Reaction<double>;

} // end namespace Antioch
