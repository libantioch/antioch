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

#ifndef ANTIOCH_CEA_THERMO_H
#define ANTIOCH_CEA_THERMO_H

// Antioch
#include "antioch_config.h"
#include "antioch/cea_curve_fit.h"
#include "antioch/cea_thermo_ascii_parsing.h"
#include "antioch/chemical_mixture.h"
#include "antioch/chemical_species.h"
#include "antioch/input_utils.h"
#include "antioch/metaprogramming.h"

// C++
#include <iomanip>
#include <vector>
#include <cmath>

namespace Antioch
{
  // Forward declarations
  template<typename CoeffType>
  class CEACurveFit;

  template<typename CoeffType=double>
  class CEAThermodynamics
  {
  public:

    CEAThermodynamics( const ChemicalMixture<CoeffType>& chem_mixture,
                       const std::string & filename = std::string(ANTIOCH_DEFAULT_FILES_PATH) + std::string(ANTIOCH_DEFAULT_THERMO_DATA) );

    //! Destructor
    /*! virtual so this can be subclassed by the user. */
    virtual ~CEAThermodynamics();

    template<typename StateType=CoeffType>
    class Cache 
    {
    public:
      const StateType &T;
      StateType T2;
      StateType T3;
      StateType T4;
      StateType lnT;
      
      explicit Cache(const StateType &T_in) 
	: T(T_in), T2(T*T), T3(T2*T), T4(T2*T2), lnT(T_in) {
        using std::log;
        lnT = log(T);
      }

      Cache(const StateType &T_in, 
            const StateType &T2_in, 
            const StateType &T3_in, 
            const StateType &T4_in, 
            const StateType &lnT_in) 
	: T(T_in), T2(T2_in), T3(T3_in), T4(T4_in), lnT(lnT_in) {
        //! \todo - correctness assertions?
      }
    private:
      Cache();      
    };
    
    void add_curve_fit( const std::string& species_name, const std::vector<CoeffType>& coeffs );

    //! Checks that curve fits have been specified for all species in the mixture.
    bool check() const;

    //! We currently need different specializations for scalar vs vector inputs here
    template<typename StateType>
    StateType
    cp( const Cache<StateType> &cache, unsigned int species ) const;

    template<typename StateType, typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value, StateType
    >::type 
    cp( const Cache<StateType> &cache, const VectorStateType& mass_fractions ) const;

    template<typename StateType>
    StateType cv( const Cache<StateType> &cache, unsigned int species ) const;

    template<typename StateType, typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value, StateType
    >::type 
    cv( const Cache<StateType> &cache, const VectorStateType& mass_fractions ) const;

    template<typename StateType>
    StateType h( const Cache<StateType> &cache, unsigned int species ) const;

    template<typename StateType, typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value, void
    >::type 
    h( const Cache<StateType> &cache, VectorStateType& h ) const;

    //! We currently need different specializations for scalar vs vector inputs here
    template<typename StateType>
    StateType
    h_RT_minus_s_R( const Cache<StateType> &cache, unsigned int species ) const;

    template<typename StateType, typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value, void
    >::type 
    h_RT_minus_s_R( const Cache<StateType> &cache, VectorStateType& h_RT_minus_s_R ) const;



    //! We currently need different specializations for scalar vs vector inputs here
    template<typename StateType>
    StateType
    dh_RT_minus_s_R_dT( const Cache<StateType> &cache, unsigned int species ) const;

    template<typename StateType, typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value, void
    >::type 
    dh_RT_minus_s_R_dT( const Cache<StateType> &cache, VectorStateType& dh_RT_minus_s_R_dT ) const;


    template<typename StateType>
    StateType cp_over_R( const Cache<StateType> &cache, unsigned int species ) const;

    template<typename StateType>
    StateType h_over_RT( const Cache<StateType> &cache, unsigned int species ) const;

    template<typename StateType>
    StateType s_over_R( const Cache<StateType> &cache, unsigned int species ) const;

    const ChemicalMixture<CoeffType>& chemical_mixture() const;

  protected:

    const ChemicalMixture<CoeffType>& _chem_mixture;

    std::vector<CEACurveFit<CoeffType>* > _species_curve_fits;

    std::vector<CoeffType> _cp_at_200p1;

  private:
    
    //! Default constructor
    /*! Private to force to user to supply a ChemicalMixture object.*/
    CEAThermodynamics();
  };

  
  /* ------------------------- Inline Functions -------------------------*/

  template<typename CoeffType>
  inline
  CEAThermodynamics<CoeffType>::CEAThermodynamics( const ChemicalMixture<CoeffType>& chem_mixture, const std::string & filename )
    : _chem_mixture(chem_mixture),
      _species_curve_fits(chem_mixture.n_species(), NULL)
  {
    // Read in CEA coefficients. Note this assumes chem_mixture is fully constructed.
    /*! \todo Generalize this to optionally read in a file instead of using the default here.
        The method is there for the reading, just need to do input file foo. */
    read_cea_thermo_data_ascii(*this, filename);

    // Cache cp values at small temperatures
    _cp_at_200p1.reserve( _species_curve_fits.size() );
    for( unsigned int s = 0; s < _species_curve_fits.size(); s++ )
      {
	_cp_at_200p1.push_back( this->cp( Cache<CoeffType>(200.1), s ) );
      }

    return;
  }


  template<typename CoeffType>
  inline
  CEAThermodynamics<CoeffType>::~CEAThermodynamics()
  {
    // Clean up all the CEACurveFits we created
    for( typename std::vector<CEACurveFit<CoeffType>* >::iterator it = _species_curve_fits.begin();
	 it < _species_curve_fits.end(); ++it )
      {
	delete (*it);
      }

    return;
  }


  template<typename CoeffType>
  inline
  void CEAThermodynamics<CoeffType>::add_curve_fit( const std::string& species_name,
						    const std::vector<CoeffType>& coeffs )
  {
    antioch_assert( _chem_mixture.species_name_map().find(species_name) !=
		    _chem_mixture.species_name_map().end() );

    unsigned int s = _chem_mixture.species_name_map().find(species_name)->second;

    antioch_assert_less_equal( s, _species_curve_fits.size() );
    antioch_assert( !_species_curve_fits[s] );

    _species_curve_fits[s] = new CEACurveFit<CoeffType>(coeffs);
    return;
  }


  template<typename CoeffType>
  inline
  bool CEAThermodynamics<CoeffType>::check() const
  {
    bool valid = true;

    for( typename std::vector<CEACurveFit<CoeffType>* >::const_iterator it = _species_curve_fits.begin();
	 it != _species_curve_fits.end(); ++ it )
      {
	if( !(*it) )
	  valid = false;
      }

    return valid;
  }


  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType
  CEAThermodynamics<CoeffType>::cp( const Cache<StateType> &cache, unsigned int species ) const
  {
    typedef typename Antioch::value_type<StateType>::type ScalarType;

    // T < 200.1 ? cp_at_200p1 : R * cp_over_R
    return
      Antioch::if_else
        (cache.T < ScalarType(200.1),
         Antioch::constant_clone
	   (cache.T,_cp_at_200p1[species]),
	 StateType
	   (this->_chem_mixture.R(species) * 
	    this->cp_over_R(cache, species)));
  }


  template<typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value, StateType
  >::type 
  CEAThermodynamics<CoeffType>::cp( const Cache<StateType> &cache,
				    const VectorStateType& mass_fractions ) const
  {
    antioch_assert_equal_to( mass_fractions.size(), _species_curve_fits.size() );
    antioch_assert_greater( mass_fractions.size(), 0 );

    StateType cp = mass_fractions[0]*this->cp(cache,0);

    for( unsigned int s = 1; s < _species_curve_fits.size(); s++ )
      {
	cp += mass_fractions[s]*this->cp(cache,s);
      }

    return cp;
  }


  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType CEAThermodynamics<CoeffType>::h( const Cache<StateType> &cache, unsigned int species ) const
  {
    return this->_chem_mixture.R(species)*cache.T*this->h_over_RT(cache,species);
  }


  template<typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value, void
  >::type 
  CEAThermodynamics<CoeffType>::h( const Cache<StateType> &cache, VectorStateType& h ) const
  {
    antioch_assert_equal_to( h.size(), _species_curve_fits.size() );

    for( unsigned int s = 0; s < _species_curve_fits.size(); s++ )
      {
	h[s] = this->_chem_mixture.R(s)*cache.T*this->h_over_RT(cache,s);
      }

    return;
  }


  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType
  CEAThermodynamics<CoeffType>::cp_over_R( const Cache<StateType> &cache, unsigned int species ) const
  {
    antioch_assert_less( species, _species_curve_fits.size() );
    // FIXME - we need assert_less to be vectorizable
    // antioch_assert_less( _species_curve_fits[species]->interval(cache.T),
    //                      _species_curve_fits[species]->n_intervals() );

    typedef typename
      Antioch::rebind<StateType, unsigned int>::type UIntType;
    const UIntType interval = this->_species_curve_fits[species]->interval(cache.T);

    const unsigned int begin_interval = Antioch::min(interval);
    const unsigned int end_interval = Antioch::max(interval)+1;
    
    // FIXME - this needs expression templates to be faster...

    StateType returnval = Antioch::zero_clone(cache.T);

    /* cp/R =  a0*T^-2   + a1*T^-1     + a2     + a3*T   + a4*T^2   + a5*T^3  + a6*T^4 */
    for (unsigned int i=begin_interval; i != end_interval; ++i)
      {
        const CoeffType * const a =
          this->_species_curve_fits[species]->coefficients(i);
	returnval = Antioch::if_else
	  (interval == i,
	   StateType(a[0]/cache.T2 + a[1]/cache.T + a[2] + a[3]*cache.T +
	             a[4]*cache.T2 + a[5]*cache.T3 + a[6]*cache.T4),
	   returnval);
      }

    return returnval;
  }


  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType CEAThermodynamics<CoeffType>::h_over_RT( const Cache<StateType> &cache, unsigned int species ) const
  {
    antioch_assert_less( species, _species_curve_fits.size() );
    antioch_assert_less( _species_curve_fits[species]->interval(cache.T),
			 _species_curve_fits[species]->n_intervals() );

    const unsigned int interval = this->_species_curve_fits[species]->interval(cache.T);
    
    const CoeffType *a = this->_species_curve_fits[species]->coefficients(interval);
    
    typedef typename
      Antioch::raw_value_type<StateType>::type raw_type;

    /* h/RT = -a0*T^-2   + a1*T^-1*lnT + a2     + a3*T/2 + a4*T^2/3 + a5*T^3/4 + a6*T^4/5 + a8/T */
    return -a[0]/cache.T2 + a[1]*cache.lnT/cache.T + a[2] +
	   a[3]*cache.T/raw_type(2) + a[4]*cache.T2/raw_type(3) +
	   a[5]*cache.T3/raw_type(4) + a[6]*cache.T4/raw_type(5) +
           a[8]/cache.T;
  }


  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType CEAThermodynamics<CoeffType>::s_over_R( const Cache<StateType> &cache, unsigned int species ) const
  {
    antioch_assert_less( species, _species_curve_fits.size() );
    antioch_assert_less( _species_curve_fits[species]->interval(cache.T),
			 _species_curve_fits[species]->n_intervals() );

    const unsigned int interval = this->_species_curve_fits[species]->interval(cache.T);
    
    const CoeffType *a = this->_species_curve_fits[species]->coefficients(interval);
    
    typedef typename
      Antioch::raw_value_type<StateType>::type raw_type;

    /* s/R = -a0*T^-2/2 - a1*T^-1     + a2*lnT + a3*T   + a4*T^2/2 + a5*T^3/3 + a6*T^4/4 + a9 */
    return -a[0]/cache.T2/raw_type(2) - a[1]/cache.T + 
	   a[2]*cache.lnT + a[3]*cache.T + a[4]*cache.T2/raw_type(2) +
           a[5]*cache.T3/raw_type(3) + a[6]*cache.T4/raw_type(4) + a[9];
  }


  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType
  CEAThermodynamics<CoeffType>::h_RT_minus_s_R( const Cache<StateType> &cache, unsigned int species ) const
  {
    antioch_assert_less( species, _species_curve_fits.size() );
    // FIXME - we need assert_less to be vectorizable
//    antioch_assert_less( _species_curve_fits[species]->interval(cache.T),
//			 _species_curve_fits[species]->n_intervals() );

    typedef typename
      Antioch::rebind<StateType, unsigned int>::type UIntType;
    const UIntType interval = this->_species_curve_fits[species]->interval(cache.T);
    const unsigned int begin_interval = Antioch::min(interval);
    const unsigned int end_interval = Antioch::max(interval)+1;
    
    // FIXME - this needs expression templates to be faster...

    StateType returnval = Antioch::zero_clone(cache.T);

    typedef typename
      Antioch::raw_value_type<StateType>::type raw_type;

    /* h/RT = -a[0]/T2    + a[1]*lnT/T + a[2]     + a[3]*T/2. + a[4]*T2/3. + a[5]*T3/4. + a[6]*T4/5. + a[8]/T,
       s/R  = -a[0]/T2/2. - a[1]/T     + a[2]*lnT + a[3]*T    + a[4]*T2/2. + a[5]*T3/3. + a[6]*T4/4. + a[9]   */
    for (unsigned int i=begin_interval; i != end_interval; ++i)
      {
        const CoeffType * const a =
          this->_species_curve_fits[species]->coefficients(i);
	returnval = Antioch::if_else
	  (interval == i,
	   StateType(-a[0]/cache.T2/raw_type(2) + 
                     (a[1] + a[8])/cache.T + 
		     a[1]*cache.lnT/cache.T - a[2]*cache.lnT + 
		     (a[2] - a[9]) - a[3]*cache.T/raw_type(2) - 
		     a[4]*cache.T2/raw_type(6) -
                     a[5]*cache.T3/raw_type(12) - 
		     a[6]*cache.T4/raw_type(20)),
	   returnval);
      }

    return returnval;
  }


  template<typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value, void
  >::type 
  CEAThermodynamics<CoeffType>::h_RT_minus_s_R( const Cache<StateType> &cache,
						VectorStateType& h_RT_minus_s_R ) const
  {
    antioch_assert_equal_to( h_RT_minus_s_R.size(), _species_curve_fits.size() );

    for( unsigned int s = 0; s < _species_curve_fits.size(); s++ )
      {
	h_RT_minus_s_R[s] = this->h_RT_minus_s_R(cache,s);
      }

    return;
  }



  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType
  CEAThermodynamics<CoeffType>::dh_RT_minus_s_R_dT( const Cache<StateType> &cache, unsigned int species ) const
  {
    antioch_assert_less( species, _species_curve_fits.size() );
    // FIXME - we need assert_less to be vectorizable
    // antioch_assert_less( _species_curve_fits[species]->interval(cache.T),
    //                      _species_curve_fits[species]->n_intervals() );
    
    typedef typename
      Antioch::rebind<StateType, unsigned int>::type UIntType;
    const UIntType interval = this->_species_curve_fits[species]->interval(cache.T);

    const unsigned int begin_interval = Antioch::min(interval);
    const unsigned int end_interval = Antioch::max(interval)+1;
    
    // FIXME - this needs expression templates to be faster...

    StateType returnval = Antioch::zero_clone(cache.T);

    typedef typename
      Antioch::raw_value_type<StateType>::type raw_type;

    /* h/RT = -a[0]/T2    + a[1]*lnT/T + a[2]     + a[3]*T/2. + a[4]*T2/3. + a[5]*T3/4. + a[6]*T4/5. + a[8]/T,
       s/R  = -a[0]/T2/2. - a[1]/T     + a[2]*lnT + a[3]*T    + a[4]*T2/2. + a[5]*T3/3. + a[6]*T4/4. + a[9]   */
    for (unsigned int i=begin_interval; i != end_interval; ++i)
      {
        const CoeffType * const a =
          this->_species_curve_fits[species]->coefficients(i);
	returnval = Antioch::if_else
	  (interval == i,
	   StateType(a[0]/cache.T3 - a[8]/cache.T2 -
		     a[1]*cache.lnT/cache.T2 - a[2]/cache.T -
		     a[3]/raw_type(2) - a[4]*cache.T/raw_type(3) -
                     a[5]*cache.T2/raw_type(4) -
		     a[6]*cache.T3/raw_type(5)),
	   returnval);
      }

    return returnval;
  }


  template<typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value, void
  >::type 
  CEAThermodynamics<CoeffType>::dh_RT_minus_s_R_dT( const Cache<StateType> &cache,
                                                    VectorStateType& dh_RT_minus_s_R_dT ) const
  {
    antioch_assert_equal_to( dh_RT_minus_s_R_dT.size(), _species_curve_fits.size() );

    for( unsigned int s = 0; s < _species_curve_fits.size(); s++ )
      {
	dh_RT_minus_s_R_dT[s] = this->dh_RT_minus_s_R_dT(cache,s);
      }

    return;
  }




  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType CEAThermodynamics<CoeffType>::cv( const Cache<StateType> &cache, unsigned int species ) const
  {
    return this->cp(cache,species) - _chem_mixture.R(species);
  }

  template<typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value, StateType
  >::type 
  CEAThermodynamics<CoeffType>::cv( const Cache<StateType> &cache,
				    const VectorStateType& mass_fractions ) const
  {
    return this->cp(cache,mass_fractions) - _chem_mixture.R(mass_fractions);
  }

  template<typename CoeffType>
  inline
  const ChemicalMixture<CoeffType>& CEAThermodynamics<CoeffType>::chemical_mixture() const
  {
    return _chem_mixture;
  }

} // end namespace Antioch

#endif //ANTIOCH_CEA_THERMO_H
