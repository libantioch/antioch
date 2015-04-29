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

#ifndef ANTIOCH_NASA_EVALUATOR_H
#define ANTIOCH_NASA_EVALUATOR_H

// Antioch
#include "antioch/chemical_mixture.h"
#include "antioch/nasa_mixture.h"
#include "antioch/temp_cache.h"
#include "antioch/thermo_enum.h"

namespace Antioch
{

  template<typename CoeffType> 
  class NASACurveFit;

  template<typename CoeffType, typename NASAFit> 
  class NASAThermoMixture;

  template<typename CoeffType=double, typename NASAFit = NASA9CurveFit<CoeffType> >
  class NASAEvaluator
  {
  public:

    NASAEvaluator( const NASAThermoMixture<CoeffType, NASAFit>& cea_mixture );
    ~NASAEvaluator();

    const NASAThermoMixture<CoeffType,NASAFit>& cea_mixture() const;

    //! We currently need different specializations for scalar vs vector inputs here
    template<typename StateType>
    StateType
    cp( const TempCache<StateType>& cache, unsigned int species ) const;

    template<typename StateType, typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value, StateType
      >::type 
    cp( const TempCache<StateType>& cache, const VectorStateType& mass_fractions ) const;

    template<typename StateType>
    ANTIOCH_AUTO(StateType)
    cv( const TempCache<StateType>& cache, unsigned int species ) const
    ANTIOCH_AUTOFUNC(StateType,
		     this->cp(cache,species) -
		     this->chem_mixture().R(species))

    template<typename StateType, typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value, StateType
      >::type 
    cv( const TempCache<StateType>& cache, const VectorStateType& mass_fractions ) const;

    template<typename StateType>
    ANTIOCH_AUTO(StateType)
    h( const TempCache<StateType>& cache, unsigned int species ) const
    ANTIOCH_AUTOFUNC(StateType,
		     this->chem_mixture().R(species) * cache.T *
		     this->h_over_RT(cache,species))

    template<typename StateType, typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value, void
      >::type 
    h( const TempCache<StateType>& cache, VectorStateType& h ) const;

    //! We currently need different specializations for scalar vs vector inputs here
    template<typename StateType>
    StateType
    h_RT_minus_s_R( const TempCache<StateType>& cache, unsigned int species ) const;

    template<typename StateType, typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value, void
      >::type 
    h_RT_minus_s_R( const TempCache<StateType>& cache, VectorStateType& h_RT_minus_s_R ) const;

    template<typename StateType>
    StateType
    dh_RT_minus_s_R_dT( const TempCache<StateType>& cache, unsigned int species ) const;

    template<typename StateType, typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value, void
      >::type 
    dh_RT_minus_s_R_dT( const TempCache<StateType>& cache, VectorStateType& h_RT_minus_s_R ) const;


     //! Cp over R, directly from the curve fit
    template<typename StateType>
    StateType cp_over_R( const TempCache<StateType>& cache, unsigned int species ) const;

     //! Cv over R, from ideal gas, 
     // \f[ \frac{C_v}{\mathrm{R} = \frac{C_p}{\mathrm{R}} - 1\f]
    template<typename StateType>
    StateType cv_over_R( const TempCache<StateType>& cache, unsigned int species ) const;

     //! h over RT, directly from the curve fit
    template<typename StateType>
    StateType h_over_RT( const TempCache<StateType>& cache, unsigned int species ) const;

     //! s over R, directly from the curve fit
    template<typename StateType>
    StateType s_over_R( const TempCache<StateType>& cache, unsigned int species ) const;

    //! sensitivity of cp for a species
    template <typename StateType>
    ANTIOCH_AUTO(StateType)
    sensitivity_cp_over_R(const TempCache<StateType> & cache, MacroThermo::Parameters parameter, unsigned int species) const;

    //! sensitivity of cv for a species
    template <typename StateType>
    ANTIOCH_AUTO(StateType)
    sensitivity_cv_over_R(const TempCache<StateType> & cache, MacroThermo::Parameters parameter, unsigned int species) const;

    //! sensitivity of h for a species
    template <typename StateType>
    ANTIOCH_AUTO(StateType)
    sensitivity_h_over_RT(const TempCache<StateType> & cache, MacroThermo::Parameters parameter, unsigned int species) const;

    //! sensitivity of s for a species
    template <typename StateType>
    ANTIOCH_AUTO(StateType)
    sensitivity_s_over_R(const TempCache<StateType> & cache, MacroThermo::Parameters parameter, unsigned int species) const;

    //! sensitivity of g for a species
    template <typename StateType>
    ANTIOCH_AUTO(StateType)
    sensitivity_h_RT_minus_s_R(const TempCache<StateType> & cache, MacroThermo::Parameters parameter, unsigned int species) const;
    

  protected:

    const NASAThermoMixture<CoeffType,NASAFit>& _nasa_mixture;

    //! Convenience function
    unsigned int n_species() const;

    //! Convenience function
    const ChemicalMixture<CoeffType>& chem_mixture() const;

  private:
    
    //! Default constructor
    /*! Private to force to user to supply a NASAThermoMixture object.*/
    NASAEvaluator();

  };

  /* --------------------- Constructor/Destructor -----------------------*/
  template<typename CoeffType, typename NASAFit>
  NASAEvaluator<CoeffType,NASAFit>::NASAEvaluator( const NASAThermoMixture<CoeffType,NASAFit>& cea_mixture )
    : _nasa_mixture(cea_mixture)
  {
    return;
  }

  template<typename CoeffType, typename NASAFit>
  NASAEvaluator<CoeffType,NASAFit>::~NASAEvaluator()
  {
    return;
  }

  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType, typename NASAFit>
  inline
  const NASAThermoMixture<CoeffType,NASAFit>& NASAEvaluator<CoeffType,NASAFit>::cea_mixture() const
  {
    return _nasa_mixture;
  }
  
  template<typename CoeffType, typename NASAFit>
  inline
  unsigned int NASAEvaluator<CoeffType,NASAFit>::n_species() const
  {
    return _nasa_mixture.chemical_mixture().n_species();
  }

  template<typename CoeffType, typename NASAFit>
  inline
  const ChemicalMixture<CoeffType>& NASAEvaluator<CoeffType,NASAFit>::chem_mixture() const
  {
    return _nasa_mixture.chemical_mixture();
  }

  template<typename CoeffType, typename NASAFit>
  template<typename StateType>
  inline
  StateType
  NASAEvaluator<CoeffType,NASAFit>::cp( const TempCache<StateType>& cache, unsigned int species ) const
  {
    typedef typename Antioch::value_type<StateType>::type ScalarType;
    // T < 200.1 ? cp_at_200p1 : R * cp_over_R
        
    return
      Antioch::if_else
        (cache.T < ScalarType(200.1),
         Antioch::constant_clone
	   (cache.T,_nasa_mixture.cp_at_200p1(species)),
	 StateType
	   (this->chem_mixture().R(species) * 
	    this->cp_over_R(cache, species)));
  }

  template<typename CoeffType, typename NASAFit>
  template<typename StateType, typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value, StateType
    >::type 
  NASAEvaluator<CoeffType,NASAFit>::cp( const TempCache<StateType>& cache,
                               const VectorStateType& mass_fractions ) const
  {
    antioch_assert_equal_to( mass_fractions.size(), this->n_species() );
    antioch_assert_greater( mass_fractions.size(), 0 );
    
    StateType cp = mass_fractions[0]*this->cp(cache,0);
    
    for( unsigned int s = 1; s < this->n_species(); s++ )
      {
        cp += mass_fractions[s]*this->cp(cache,s);
      }
    
    return cp;
  }


  template<typename CoeffType, typename NASAFit>
  template<typename StateType, typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value, void
    >::type 
  NASAEvaluator<CoeffType,NASAFit>::h( const TempCache<StateType>& cache, VectorStateType& h ) const
  {
    antioch_assert_equal_to( h.size(), this->n_species() );
    
    for( unsigned int s = 0; s < this->n_species(); s++ )
      {
        h[s] = this->chem_mixture().R(s)*cache.T*this->h_over_RT(cache,s);
      }
    
    return;
  }


  template<typename CoeffType, typename NASAFit>
  template<typename StateType>
  inline
  StateType
  NASAEvaluator<CoeffType,NASAFit>::cp_over_R( const TempCache<StateType>& cache, unsigned int species ) const
  {
    antioch_assert_less( species, this->n_species() );
    // FIXME - we need assert_less to be vectorizable
    // antioch_assert_less( _nasa_mixture.curve_fit(species).interval(cache.T),
    //                      _nasa_mixture.curve_fit(species).n_intervals() );

    return this->_nasa_mixture.curve_fit(species).cp_over_R(cache);
  }

  template<typename CoeffType, typename NASAFit>
  template<typename StateType>
  inline
  StateType
  NASAEvaluator<CoeffType,NASAFit>::cv_over_R( const TempCache<StateType>& cache, unsigned int species ) const
  {
    antioch_assert_less( species, this->n_species() );
    // FIXME - we need assert_less to be vectorizable
    // antioch_assert_less( _nasa_mixture.curve_fit(species).interval(cache.T),
    //                      _nasa_mixture.curve_fit(species).n_intervals() );
    return this->_nasa_mixture.curve_fit(species).cp_over_R(cache) - StateType(1.L);
  }

  template<typename CoeffType, typename NASAFit>
  template<typename StateType>
  inline
  StateType NASAEvaluator<CoeffType,NASAFit>::h_over_RT( const TempCache<StateType>& cache, unsigned int species ) const
  {
    antioch_assert_less( species, this->n_species() );
    antioch_assert_less( _nasa_mixture.curve_fit(species).interval(cache.T),
                         _nasa_mixture.curve_fit(species).n_intervals() );
    
    return this->_nasa_mixture.curve_fit(species).h_over_RT(cache);
  }


  template<typename CoeffType, typename NASAFit>
  template<typename StateType>
  inline
  StateType NASAEvaluator<CoeffType,NASAFit>::s_over_R( const TempCache<StateType>& cache, unsigned int species ) const
  {
    antioch_assert_less( species, this->n_species() );
    antioch_assert_less( _nasa_mixture.curve_fit(species).interval(cache.T),
                         _nasa_mixture.curve_fit(species).n_intervals() );
    
    return this->_nasa_mixture.curve_fit(species).s_over_R(cache);
  }
  

  template<typename CoeffType, typename NASAFit>
  template<typename StateType>
  inline
  StateType
  NASAEvaluator<CoeffType,NASAFit>::h_RT_minus_s_R( const TempCache<StateType>& cache, unsigned int species ) const
  {
    antioch_assert_less( species, this->n_species() );
    // FIXME - we need assert_less to be vectorizable
    // antioch_assert_less( _nasa_mixture.curve_fit(species).interval(cache.T),
    //                      _nasa_mixture.curve_fit(species).n_intervals() );
    
    return this->_nasa_mixture.curve_fit(species).h_RT_minus_s_R(cache);
  }


  template<typename CoeffType, typename NASAFit>
  template<typename StateType, typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value, void
    >::type 
  NASAEvaluator<CoeffType,NASAFit>::h_RT_minus_s_R( const TempCache<StateType>& cache,
                                           VectorStateType& h_RT_minus_s_R ) const
  {
    antioch_assert_equal_to( h_RT_minus_s_R.size(), this->n_species() );
    
    for( unsigned int s = 0; s < this->n_species(); s++ )
      {
        h_RT_minus_s_R[s] = this->h_RT_minus_s_R(cache,s);
      }
    
    return;
  }



  template<typename CoeffType, typename NASAFit>
  template<typename StateType>
  inline
  StateType
  NASAEvaluator<CoeffType,NASAFit>::dh_RT_minus_s_R_dT( const TempCache<StateType>& cache, unsigned int species ) const
  {
    antioch_assert_less( species, this->n_species() );
    // FIXME - we need assert_less to be vectorizable
    // antioch_assert_less( _nasa_mixture.curve_fit(species).interval(cache.T),
    //                      _nasa_mixture.curve_fit(species).n_intervals() );
      
    return this->_nasa_mixture.curve_fit(species).dh_RT_minus_s_R_dT(cache);
  }


  template<typename CoeffType, typename NASAFit>
  template<typename StateType, typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value, void
    >::type 
  NASAEvaluator<CoeffType,NASAFit>::dh_RT_minus_s_R_dT( const TempCache<StateType>& cache,
                                               VectorStateType& dh_RT_minus_s_R_dT ) const
  {
    antioch_assert_equal_to( dh_RT_minus_s_R_dT.size(), this->n_species() );
    
    for( unsigned int s = 0; s < this->n_species(); s++ )
      {
        dh_RT_minus_s_R_dT[s] = this->dh_RT_minus_s_R_dT(cache,s);
      }
    
    return;
  }
  
  
  template<typename CoeffType, typename NASAFit>
  template<typename StateType, typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value, StateType
    >::type 
  NASAEvaluator<CoeffType,NASAFit>::cv( const TempCache<StateType>& cache,
                               const VectorStateType& mass_fractions ) const
    {
      return this->cp(cache,mass_fractions) - this->chem_mixture().R(mass_fractions);
    }

  template<typename CoeffType, typename NASAFit>
  template <typename StateType>
  inline
  ANTIOCH_AUTO(StateType)
  NASAEvaluator<CoeffType,NASAFit>::sensitivity_cp_over_R(const TempCache<StateType> & cache, MacroThermo::Parameters parameter, unsigned int species) const
  {
      antioch_assert_less(species,this->n_species());
      return _nasa_mixture.curve_fit(species).dcp_over_R_dpar(cache,parameter);
  }

  template<typename CoeffType, typename NASAFit>
  template <typename StateType>
  inline
  ANTIOCH_AUTO(StateType)
  NASAEvaluator<CoeffType,NASAFit>::sensitivity_cv_over_R(const TempCache<StateType> & cache, MacroThermo::Parameters parameter, unsigned int species) const
  {
      antioch_assert_less(species,this->n_species());
      return _nasa_mixture.curve_fit(species).dcp_over_R_dpar(cache,parameter);
  }

  template<typename CoeffType, typename NASAFit>
  template <typename StateType>
  inline
  ANTIOCH_AUTO(StateType)
  NASAEvaluator<CoeffType,NASAFit>::sensitivity_h_over_RT(const TempCache<StateType> & cache, MacroThermo::Parameters parameter, unsigned int species) const
  {
      antioch_assert_less(species,this->n_species());
      return _nasa_mixture.curve_fit(species).dh_over_RT_dpar(cache,parameter);
  }

  template<typename CoeffType, typename NASAFit>
  template <typename StateType>
  inline
  ANTIOCH_AUTO(StateType)
  NASAEvaluator<CoeffType,NASAFit>::sensitivity_s_over_R(const TempCache<StateType> & cache, MacroThermo::Parameters parameter, unsigned int species) const
  {
      antioch_assert_less(species,this->n_species());
      return _nasa_mixture.curve_fit(species).ds_over_R_dpar(cache,parameter);
  }

  template<typename CoeffType, typename NASAFit>
  template <typename StateType>
  inline
  ANTIOCH_AUTO(StateType)
  NASAEvaluator<CoeffType,NASAFit>::sensitivity_h_RT_minus_s_R(const TempCache<StateType> & cache, MacroThermo::Parameters parameter, unsigned int species) const
  {
      antioch_assert_less(species,this->n_species());
      return _nasa_mixture.curve_fit(species).dh_RT_minus_s_R_dpar(cache,parameter);
  }

} // end namespace Antioch

#endif // ANTIOCH_NASA_EVALUATOR_H
