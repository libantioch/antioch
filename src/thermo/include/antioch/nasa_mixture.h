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

#ifndef ANTIOCH_NASA_MIXTURE_H
#define ANTIOCH_NASA_MIXTURE_H

// Antioch
#include "antioch/chemical_mixture.h"
#include "antioch/chemical_species.h"
#include "antioch/input_utils.h"
#include "antioch/cea_curve_fit.h"
#include "antioch/nasa_curve_fit.h"
#include "antioch/temp_cache.h"

// C++
#include <iomanip>
#include <vector>
#include <cmath>

namespace Antioch
{
  
  template<typename CoeffType, typename NASAFit>
  class NASAEvaluator;

  template<typename CoeffType=double, typename NASAFit = NASA9CurveFit<CoeffType> >
  class NASAThermoMixture
  {
  public:

    NASAThermoMixture( const ChemicalMixture<CoeffType>& chem_mixture );

    //! Destructor
    /*! virtual so this can be subclassed by the user. */
    virtual ~NASAThermoMixture();

     //NASA9 friendly
    void add_curve_fit( const std::string& species_name, const std::vector<CoeffType>& coeffs );

     //NASA general
    void add_curve_fit( const std::string& species_name, const std::vector<CoeffType>& coeffs, const std::vector<CoeffType>& temps );

    const NASAFit& curve_fit( unsigned int s ) const;

    CoeffType cp_at_200p1( unsigned int s ) const;

    //! Checks that curve fits have been specified for all species in the mixture.
    bool check() const;

    const ChemicalMixture<CoeffType>& chemical_mixture() const;
    
  protected:

    const ChemicalMixture<CoeffType>& _chem_mixture;

    std::vector<NASAFit* > _species_curve_fits;

    std::vector<CoeffType> _cp_at_200p1;

  private:
    
    //! Default constructor
    /*! Private to force to user to supply a ChemicalMixture object.*/
    NASAThermoMixture();

  };

  /* --------------------- Constructor/Destructor -----------------------*/
  template<typename CoeffType, typename NASAFit>
  NASAThermoMixture<CoeffType,NASAFit>::NASAThermoMixture( const ChemicalMixture<CoeffType>& chem_mixture )
    : _chem_mixture(chem_mixture),
      _species_curve_fits(chem_mixture.n_species(), NULL),
      _cp_at_200p1( _species_curve_fits.size() )
  {
    return;
  }


  template<typename CoeffType, typename NASAFit>
  NASAThermoMixture<CoeffType,NASAFit>::~NASAThermoMixture()
  {
    // Clean up all the NASAFits we created
    for( typename std::vector<NASAFit* >::iterator it = _species_curve_fits.begin();
	 it < _species_curve_fits.end(); ++it )
      {
	delete (*it);
      }

    return;
  }

  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType, typename NASAFit>
  inline
  void NASAThermoMixture<CoeffType,NASAFit>::add_curve_fit( const std::string& species_name,
                                                   const std::vector<CoeffType>& coeffs)
  {
    antioch_assert( _chem_mixture.species_name_map().find(species_name) !=
		    _chem_mixture.species_name_map().end() );

    unsigned int s = _chem_mixture.species_name_map().find(species_name)->second;

    antioch_assert_less_equal( s, _species_curve_fits.size() );
    antioch_assert( !_species_curve_fits[s] );

    _species_curve_fits[s] = new NASAFit(coeffs);

    NASAEvaluator<CoeffType,NASAFit> evaluator( *this );
    _cp_at_200p1[s] = evaluator.cp( TempCache<CoeffType>(200.1), s );

    return;
  }

  template<typename CoeffType, typename NASAFit>
  inline
  void NASAThermoMixture<CoeffType,NASAFit>::add_curve_fit( const std::string& species_name,
                                                   const std::vector<CoeffType>& coeffs, const std::vector<CoeffType>& temps )
  {
    antioch_assert( _chem_mixture.species_name_map().find(species_name) !=
		    _chem_mixture.species_name_map().end() );

    unsigned int s = _chem_mixture.species_name_map().find(species_name)->second;

    antioch_assert_less_equal( s, _species_curve_fits.size() );
    antioch_assert( !_species_curve_fits[s] );

    _species_curve_fits[s] = new NASAFit(coeffs,temps);

    NASAEvaluator<CoeffType,NASAFit> evaluator( *this );
    _cp_at_200p1[s] = evaluator.cp( TempCache<CoeffType>(200.1), s );

    return;
  }


  template<typename CoeffType, typename NASAFit>
  inline
  bool NASAThermoMixture<CoeffType,NASAFit>::check() const
  {
    bool valid = true;

    for( typename std::vector<NASAFit* >::const_iterator it = _species_curve_fits.begin();
	 it != _species_curve_fits.end(); ++ it )
      {
	if( !(*it) )
	  valid = false;
      }

    return valid;
  }

  template<typename CoeffType, typename NASAFit>
  inline
  const NASAFit& NASAThermoMixture<CoeffType,NASAFit>::curve_fit( unsigned int s ) const
  {
    antioch_assert_less( s, _species_curve_fits.size() );
    return *_species_curve_fits[s];
  }

  template<typename CoeffType, typename NASAFit>
  inline
  CoeffType NASAThermoMixture<CoeffType,NASAFit>::cp_at_200p1( unsigned int s ) const
  {
    antioch_assert_less( s, _cp_at_200p1.size() );
    return _cp_at_200p1[s];
  }

  template<typename CoeffType, typename NASAFit>
  inline
  const ChemicalMixture<CoeffType>& NASAThermoMixture<CoeffType,NASAFit>::chemical_mixture() const
  {
    return _chem_mixture;
  }

} // end namespace Antioch

#endif // ANTIOCH_NASA_MIXTURE_H
