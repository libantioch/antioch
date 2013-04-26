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

#ifndef ANTIOCH_CEA_EVALUATOR_H
#define ANTIOCH_CEA_EVALUATOR_H

// Antioch
#include "antioch/chemical_mixture.h"
#include "antioch/cea_mixture.h"
#include "antioch/temp_cache.h"

namespace Antioch
{

  template<typename CoeffType> class CEAThermoMixture;

  template<typename CoeffType=double>
  class CEAEvaluator
  {
  public:

    CEAEvaluator( const CEAThermoMixture<CoeffType>& cea_mixture );
    ~CEAEvaluator();

    const CEAThermoMixture<CoeffType>& cea_mixture() const;

    //! We currently need different specializations for scalar vs vector inputs here
    template<typename StateType>
    typename enable_if_c<
      !has_size<StateType>::value, StateType
      >::type 
    cp( const TempCache<StateType>& cache, unsigned int species ) const;

    template<typename StateType>
    typename enable_if_c<
      has_size<StateType>::value, StateType
      >::type 
    cp( const TempCache<StateType>& cache, unsigned int species ) const;

    template<typename StateType, typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value, StateType
      >::type 
    cp( const TempCache<StateType>& cache, const VectorStateType& mass_fractions ) const;

    template<typename StateType>
    StateType cv( const TempCache<StateType>& cache, unsigned int species ) const;

    template<typename StateType, typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value, StateType
      >::type 
    cv( const TempCache<StateType>& cache, const VectorStateType& mass_fractions ) const;

    template<typename StateType>
    StateType h( const TempCache<StateType>& cache, unsigned int species ) const;

    template<typename StateType, typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value, void
      >::type 
    h( const TempCache<StateType>& cache, VectorStateType& h ) const;

    //! We currently need different specializations for scalar vs vector inputs here
    template<typename StateType>
    typename enable_if_c<
      !has_size<StateType>::value, StateType
      >::type 
    h_RT_minus_s_R( const TempCache<StateType>& cache, unsigned int species ) const;

    template<typename StateType>
    typename enable_if_c<
      has_size<StateType>::value, StateType
      >::type 
    h_RT_minus_s_R( const TempCache<StateType>& cache, unsigned int species ) const;

    template<typename StateType, typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value, void
      >::type 
    h_RT_minus_s_R( const TempCache<StateType>& cache, VectorStateType& h_RT_minus_s_R ) const;

    template<typename StateType>
    StateType cp_over_R( const TempCache<StateType>& cache, unsigned int species ) const;

    template<typename StateType>
    StateType h_over_RT( const TempCache<StateType>& cache, unsigned int species ) const;

    template<typename StateType>
    StateType s_over_R( const TempCache<StateType>& cache, unsigned int species ) const;

  protected:

    const CEAThermoMixture<CoeffType>& _cea_mixture;

    //! Convenience function
    unsigned int n_species() const;

    //! Convenience function
    const ChemicalMixture<CoeffType>& chem_mixture() const;

  private:
    
    //! Default constructor
    /*! Private to force to user to supply a CEAThermoMixture object.*/
    CEAEvaluator();

  };

  /* --------------------- Constructor/Destructor -----------------------*/
  template<typename CoeffType>
  CEAEvaluator<CoeffType>::CEAEvaluator( const CEAThermoMixture<CoeffType>& cea_mixture )
    : _cea_mixture(cea_mixture)
  {
    return;
  }

  template<typename CoeffType>
  CEAEvaluator<CoeffType>::~CEAEvaluator()
  {
    return;
  }

  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType>
  inline
  const CEAThermoMixture<CoeffType>& CEAEvaluator<CoeffType>::cea_mixture() const
  {
    return _cea_mixture;
  }
  
  template<typename CoeffType>
  inline
  unsigned int CEAEvaluator<CoeffType>::n_species() const
  {
    return _cea_mixture.chemical_mixture().n_species();
  }

  template<typename CoeffType>
  inline
  const ChemicalMixture<CoeffType>& CEAEvaluator<CoeffType>::chem_mixture() const
  {
    return _cea_mixture.chemical_mixture();
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  typename enable_if_c<
    !has_size<StateType>::value, StateType
    >::type 
    CEAEvaluator<CoeffType>::cp( const TempCache<StateType>& cache, unsigned int species ) const
    {
      StateType cp = 0.0;

      if( cache.T < 200.1 )
        {
          cp =  _cea_mixture.cp_at_200p1(species);
        }
      else
        {
          cp = this->chem_mixture().R(species)*this->cp_over_R(cache, species);
        }
    
      return cp;
    }

  template<typename CoeffType>
  template<typename StateType>
  inline
  typename enable_if_c<
    has_size<StateType>::value, StateType
    >::type 
  CEAEvaluator<CoeffType>::cp( const TempCache<StateType>& cache, unsigned int species ) const
  {
    // Use an input datum to make sure we get the size right
    StateType cp = Antioch::zero_clone(cache.T);
    
    const std::size_t size = cache.T.size();
    for (std::size_t i = 0; i != size; ++i)
      {
        typedef TempCache<typename Antioch::value_type<StateType>::type> SubCache;
        cp[i] = this->cp
          (SubCache(cache.T[i],cache.T2[i],cache.T3[i],cache.T4[i],cache.lnT[i]),
           species);
      }
    
    return cp;
  }
  
  template<typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value, StateType
    >::type 
  CEAEvaluator<CoeffType>::cp( const TempCache<StateType>& cache,
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


  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType CEAEvaluator<CoeffType>::h( const TempCache<StateType>& cache, unsigned int species ) const
  {
    return this->_chem_mixture.R(species)*cache.T*this->h_over_RT(cache,species);
  }


  template<typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value, void
    >::type 
  CEAEvaluator<CoeffType>::h( const TempCache<StateType>& cache, VectorStateType& h ) const
  {
    antioch_assert_equal_to( h.size(), this->n_species() );
    
    for( unsigned int s = 0; s < this->n_species(); s++ )
      {
        h[s] = this->_chem_mixture.R(s)*cache.T*this->h_over_RT(cache,s);
      }
    
    return;
  }


  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType CEAEvaluator<CoeffType>::cp_over_R( const TempCache<StateType>& cache, unsigned int species ) const
  {
    antioch_assert_less( species, this->n_species() );
    antioch_assert_less( _cea_mixture.curve_fit(species).interval(cache.T),
                         _cea_mixture.curve_fit(species).n_intervals() );
    
    const unsigned int interval = this->_cea_mixture.curve_fit(species).interval(cache.T);
    
    const CoeffType *a = this->_cea_mixture.curve_fit(species).coefficients(interval);
    
    /* cp/R =  a0*T^-2   + a1*T^-1     + a2     + a3*T   + a4*T^2   + a5*T^3  + a6*T^4 */
    return a[0]/cache.T2 + a[1]/cache.T + a[2] + a[3]*cache.T + a[4]*cache.T2 + a[5]*cache.T3 + a[6]*cache.T4;
  }


  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType CEAEvaluator<CoeffType>::h_over_RT( const TempCache<StateType>& cache, unsigned int species ) const
  {
    antioch_assert_less( species, this->n_species() );
    antioch_assert_less( _cea_mixture.curve_fit(species).interval(cache.T),
                         _cea_mixture.curve_fit(species).n_intervals() );
    
    const unsigned int interval = this->_cea_mixture.curve_fit(species).interval(cache.T);
    
    const CoeffType *a = this->_cea_mixture.curve_fit(species).coefficients(interval);
    
    /* h/RT = -a0*T^-2   + a1*T^-1*lnT + a2     + a3*T/2 + a4*T^2/3 + a5*T^3/4 + a6*T^4/5 + a8/T */
    return -a[0]/cache.T2 + a[1]*cache.lnT/cache.T + a[2] + a[3]*cache.T/2.0 + a[4]*cache.T2/3.0 + a[5]*cache.T3/4.0 + a[6]*cache.T4/5.0 + a[8]/cache.T;
  }


  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType CEAEvaluator<CoeffType>::s_over_R( const TempCache<StateType>& cache, unsigned int species ) const
  {
    antioch_assert_less( species, this->n_species() );
    antioch_assert_less( _cea_mixture.curve_fit(species).interval(cache.T),
                         _cea_mixture.curve_fit(species).n_intervals() );
    
    const unsigned int interval = this->_cea_mixture.curve_fit(species).interval(cache.T);
    
    const CoeffType *a = this->_cea_mixture.curve_fit(species).coefficients(interval);
    
    /* s/R = -a0*T^-2/2 - a1*T^-1     + a2*lnT + a3*T   + a4*T^2/2 + a5*T^3/3 + a6*T^4/4 + a9 */
    return -a[0]/cache.T2/2.0 - a[1]/cache.T + a[2]*cache.lnT + a[3]*cache.T + a[4]*cache.T2/2.0 + a[5]*cache.T3/3.0 + a[6]*cache.T4/4.0 + a[9];
  }
  

  template<typename CoeffType>
  template<typename StateType>
  inline
  typename enable_if_c<
    !has_size<StateType>::value, StateType
    >::type 
    CEAEvaluator<CoeffType>::h_RT_minus_s_R( const TempCache<StateType>& cache, unsigned int species ) const
    {
      antioch_assert_less( species, this->n_species() );
      antioch_assert_less( _cea_mixture.curve_fit(species).interval(cache.T),
                           _cea_mixture.curve_fit(species).n_intervals() );
      
      const unsigned int interval = this->_cea_mixture.curve_fit(species).interval(cache.T);
      
      const CoeffType *a = this->_cea_mixture.curve_fit(species).coefficients(interval);
      
      /* h/RT = -a[0]/T2    + a[1]*lnT/T + a[2]     + a[3]*T/2. + a[4]*T2/3. + a[5]*T3/4. + a[6]*T4/5. + a[8]/T,
         s/R  = -a[0]/T2/2. - a[1]/T     + a[2]*lnT + a[3]*T    + a[4]*T2/2. + a[5]*T3/3. + a[6]*T4/4. + a[9]   */
      return -a[0]/cache.T2/2.0 + (a[1] + a[8])/cache.T + a[1]*cache.lnT/cache.T - a[2]*cache.lnT + (a[2] - a[9]) - a[3]*cache.T/2.0 - a[4]*cache.T2/6.0 - a[5]*cache.T3/12.0 - a[6]*cache.T4/20.0;
    }


  template<typename CoeffType>
  template<typename StateType>
  inline
  typename enable_if_c<
    has_size<StateType>::value, StateType
    >::type 
  CEAEvaluator<CoeffType>::h_RT_minus_s_R( const TempCache<StateType>& cache, unsigned int species ) const
  {
    antioch_assert_less( species, this->n_species() );
    
    const std::size_t size = cache.T.size();
    
    // Use an input variable to determine sizing
    StateType returnval = Antioch::zero_clone(cache.T);
    
    for (std::size_t i = 0; i != size; ++i)
      {
        typedef TempCache<typename Antioch::value_type<StateType>::type> SubCache;
        returnval[i] = this->h_RT_minus_s_R
          (SubCache(cache.T[i],cache.T2[i],cache.T3[i],cache.T4[i],cache.lnT[i]),
           species);
      }
    
    return returnval;
  }
  
  
  template<typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value, void
    >::type 
  CEAEvaluator<CoeffType>::h_RT_minus_s_R( const TempCache<StateType>& cache,
                                           VectorStateType& h_RT_minus_s_R ) const
  {
    antioch_assert_equal_to( h_RT_minus_s_R.size(), this->n_species() );
    
    for( unsigned int s = 0; s < this->n_species(); s++ )
      {
        h_RT_minus_s_R[s] = this->h_RT_minus_s_R(cache,s);
      }
    
    return;
  }
  
  
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType CEAEvaluator<CoeffType>::cv( const TempCache<StateType>& cache, unsigned int species ) const
  {
    return this->cp(cache,species) - this->chem_mixture().R(species);
  }
  
  template<typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value, StateType
    >::type 
  CEAEvaluator<CoeffType>::cv( const TempCache<StateType>& cache,
                               const VectorStateType& mass_fractions ) const
    {
      return this->cp(cache,mass_fractions) - this->chem_mixture().R(mass_fractions);
    }

} // end namespace Antioch

#endif // ANTIOCH_CEA_EVALUATOR_H
