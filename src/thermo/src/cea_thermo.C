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

// C++
#include <cmath>
#include <sstream>
#include <limits>

// This class
#include "antioch/cea_thermo.h"

// Antioch
#include "antioch/input_utils.h"
#include "antioch/chemical_species.h"
#include "antioch/chemical_mixture.h"
#include "antioch/cea_curve_fit.h"
#include "antioch/cea_thermo_ascii_parsing.h"

namespace Antioch
{
  template<class NumericType>
  CEAThermodynamics<NumericType>::CEAThermodynamics( const ChemicalMixture<NumericType>& chem_mixture )
    : _chem_mixture(chem_mixture),
      _species_curve_fits(chem_mixture.n_species(), NULL),
      _T(std::numeric_limits<NumericType>::max())
  {
    // Read in CEA coefficients. Note this assumes chem_mixture is fully constructed.
    /*! \todo Generalize this to optionally read in a file instead of using the default here.
        The method is there for the reading, just need to do input file foo. */
    read_cea_thermo_data_ascii_default(*this);

    // Cache cp values at small temperatures
    _cp_at_200p1.reserve( _species_curve_fits.size() );
    for( unsigned int s = 0; s < _species_curve_fits.size(); s++ )
      {
	_cp_at_200p1.push_back( this->cp( 200.1, s ) );
      }

    return;
  }

  template<class NumericType>
  CEAThermodynamics<NumericType>::~CEAThermodynamics()
  {
    // Clean up all the CEACurveFits we created
    for( typename std::vector<CEACurveFit<NumericType>* >::iterator it = _species_curve_fits.begin();
	 it < _species_curve_fits.end(); ++it )
      {
	delete (*it);
      }

    return;
  }

  template<class NumericType>
  void CEAThermodynamics<NumericType>::add_curve_fit( const std::string& species_name,
						      const std::vector<NumericType>& coeffs )
  {
    antioch_assert( _chem_mixture.active_species_name_map().find(species_name) !=
		    _chem_mixture.active_species_name_map().end() );

    unsigned int s = _chem_mixture.active_species_name_map().find(species_name)->second;

    antioch_assert_less_equal( s, _species_curve_fits.size() );
    antioch_assert( !_species_curve_fits[s] );

    _species_curve_fits[s] = new CEACurveFit<NumericType>(coeffs);
    return;
  }

  template<class NumericType>
  bool CEAThermodynamics<NumericType>::check() const
  {
    bool valid = true;

    for( typename std::vector<CEACurveFit<NumericType>* >::const_iterator it = _species_curve_fits.begin();
	 it != _species_curve_fits.end(); ++ it )
      {
	if( !(*it) )
	  valid = false;
      }

    return valid;
  }

  template<class NumericType>
  NumericType CEAThermodynamics<NumericType>::cp( NumericType T, unsigned int species ) const
  {
    NumericType cp = 0.0;

    if( T < 200.1 )
      {
	cp =  _cp_at_200p1[species];
      }
    else
      {
	cp = this->_chem_mixture.R(species)*this->cp_over_R(T,species);
      }
    
    return cp;
  }

  template<class NumericType>
  NumericType CEAThermodynamics<NumericType>::cp( NumericType T,
						  const std::vector<NumericType>& mass_fractions ) const
  {
    antioch_assert_equal_to( mass_fractions.size(), _species_curve_fits.size() );

    NumericType cp = 0.0;

    for( unsigned int s = 0; s < _species_curve_fits.size(); s++ )
      {
	cp += mass_fractions[s]*this->cp(T,s);
      }

    return cp;
  }

  template<class NumericType>
  NumericType CEAThermodynamics<NumericType>::h( NumericType T, unsigned int species ) const
  {
    return this->_chem_mixture.R(species)*T*this->h_over_RT(T,species);
  }

  template<class NumericType>
  void CEAThermodynamics<NumericType>::h( NumericType T, std::vector<NumericType>& h ) const
  {
    antioch_assert_equal_to( h.size(), _species_curve_fits.size() );

    for( unsigned int s = 0; s < _species_curve_fits.size(); s++ )
      {
	h[s] = this->_chem_mixture.R(s)*T*this->h_over_RT(T,s);
      }

    return;
  }

  template<class NumericType>
  NumericType CEAThermodynamics<NumericType>::cp_over_R( NumericType T, unsigned int species ) const
  {
    antioch_assert_less( species, _species_curve_fits.size() );
    antioch_assert_less( _species_curve_fits[species]->interval(T),
			 _species_curve_fits[species]->n_intervals() );

    const unsigned int interval = this->_species_curve_fits[species]->interval(T);
    
    const NumericType *a = this->_species_curve_fits[species]->coefficients(interval);
    
    this->update_cache(T);

    /* cp/R =  a0*T^-2   + a1*T^-1     + a2     + a3*T   + a4*T^2   + a5*T^3  + a6*T^4 */
    return a[0]/_T2 + a[1]/_T + a[2] + a[3]*_T + a[4]*_T2 + a[5]*_T3 + a[6]*_T4;
  }

  template<class NumericType>
  NumericType CEAThermodynamics<NumericType>::h_over_RT( NumericType T, unsigned int species ) const
  {
    antioch_assert_less( species, _species_curve_fits.size() );
    antioch_assert_less( _species_curve_fits[species]->interval(T),
			 _species_curve_fits[species]->n_intervals() );

    const unsigned int interval = this->_species_curve_fits[species]->interval(T);
    
    const NumericType *a = this->_species_curve_fits[species]->coefficients(interval);
    
    this->update_cache(T);

    /* h/RT = -a0*T^-2   + a1*T^-1*lnT + a2     + a3*T/2 + a4*T^2/3 + a5*T^3/4 + a6*T^4/5 + a8/T */
    return -a[0]/_T2 + a[1]*_lnT/_T + a[2] + a[3]*_T/2.0 + a[4]*_T2/3.0 + a[5]*_T3/4.0 + a[6]*_T4/5.0 + a[8]/_T;
  }

  template<class NumericType>
  NumericType CEAThermodynamics<NumericType>::s_over_R( NumericType T, unsigned int species ) const
  {
    antioch_assert_less( species, _species_curve_fits.size() );
    antioch_assert_less( _species_curve_fits[species]->interval(T),
			 _species_curve_fits[species]->n_intervals() );

    const unsigned int interval = this->_species_curve_fits[species]->interval(T);
    
    const NumericType *a = this->_species_curve_fits[species]->coefficients(interval);
    
    this->update_cache(T);

    /* s/R = -a0*T^-2/2 - a1*T^-1     + a2*lnT + a3*T   + a4*T^2/2 + a5*T^3/3 + a6*T^4/4 + a9 */
    return -a[0]/_T2/2.0 - a[1]/_T + a[2]*_lnT + a[3]*_T + a[4]*_T2/2.0 + a[5]*_T3/3.0 + a[6]*_T4/4.0 + a[9];
  }

  template<class NumericType>
  NumericType CEAThermodynamics<NumericType>::h_RT_minus_s_R( NumericType T, unsigned int species ) const
  {
    antioch_assert_less( species, _species_curve_fits.size() );
    antioch_assert_less( _species_curve_fits[species]->interval(T),
			 _species_curve_fits[species]->n_intervals() );

    const unsigned int interval = this->_species_curve_fits[species]->interval(T);
    
    const NumericType *a = this->_species_curve_fits[species]->coefficients(interval);
    
    this->update_cache(T);

    /* h/RT = -a[0]/T2    + a[1]*lnT/T + a[2]     + a[3]*T/2. + a[4]*T2/3. + a[5]*T3/4. + a[6]*T4/5. + a[8]/T,
       s/R  = -a[0]/T2/2. - a[1]/T     + a[2]*lnT + a[3]*T    + a[4]*T2/2. + a[5]*T3/3. + a[6]*T4/4. + a[9]   */
    return -a[0]/_T2/2.0 + (a[1] + a[8])/_T + a[1]*_lnT/_T - a[2]*_lnT + (a[2] - a[9]) - a[3]*_T/2.0 - a[4]*_T2/6.0 - a[5]*_T3/12.0 - a[6]*_T4/20.0;
  }

  template<class NumericType>
  void CEAThermodynamics<NumericType>::h_RT_minus_s_R( NumericType T,
						       std::vector<NumericType>& h_RT_minus_s_R ) const
  {
    antioch_assert_equal_to( h_RT_minus_s_R.size(), _species_curve_fits.size() );

    for( unsigned int s = 0; s < _species_curve_fits.size(); s++ )
      {
	h_RT_minus_s_R[s] = this->h_RT_minus_s_R(T,s);
      }

    return;
  }

  /* ------------------------- Instantiate ------------------------- */
  template class CEAThermodynamics<double>;

} // end namespace Antioch
