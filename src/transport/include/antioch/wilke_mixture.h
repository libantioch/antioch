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

#ifndef ANTIOCH_WILKE_MIXTURE_H
#define ANTIOCH_WILKE_MIXTURE_H

// C++
#include <vector>

namespace Antioch
{
  // Forward declarations
  template<class CoeffType>
  class ChemicalMixture;

  template<class CoeffType=double>
  class WilkeMixture
  {
  public:

    WilkeMixture( const ChemicalMixture<CoeffType>& chem_mixture );
    ~WilkeMixture();

    CoeffType Mr_Ms_to_the_one_fourth( const unsigned int r,
                                       const unsigned int s ) const;
    
    CoeffType denominator( const unsigned int r,
                           const unsigned int s ) const;

    const ChemicalMixture<CoeffType>& chem_mixture() const;

  protected:

    const ChemicalMixture<CoeffType>& _chem_mixture;

    //! Cache for numerator term
    /*! \todo We should use a more efficient data structure */
    std::vector<std::vector<CoeffType> > _Mr_Ms_to_the_one_fourth;
    
    //! Cache for denominator term
    /*! \todo We should use a more efficient data structure */
    std::vector<std::vector<CoeffType> > _denom;

  };

  template<class CoeffType>
  WilkeMixture<CoeffType>::WilkeMixture( const ChemicalMixture<CoeffType>& chem_mixture )
    : _chem_mixture(chem_mixture),
      _Mr_Ms_to_the_one_fourth(chem_mixture.n_species()),
      _denom(chem_mixture.n_species())
  {

    for( unsigned int r = 0; r < chem_mixture.n_species(); r++ )
      {
        _Mr_Ms_to_the_one_fourth[r].resize(chem_mixture.n_species());
        _denom[r].resize(chem_mixture.n_species());

        for( unsigned int s = 0; s < chem_mixture.n_species(); s++ )
          {
            const CoeffType Mr = chem_mixture.M(r);
            const CoeffType Ms = chem_mixture.M(s);

            _Mr_Ms_to_the_one_fourth[r][s] = std::pow( Mr/Ms, 0.25 );
            _denom[r][s] = std::sqrt(8.0*(1.0+Ms/Mr));
          }
      }

    return;
  }

  template<class CoeffType>
  WilkeMixture<CoeffType>::~WilkeMixture()
  {
    return;
  }

  
  template<class CoeffType>
  inline
  CoeffType WilkeMixture<CoeffType>::Mr_Ms_to_the_one_fourth( const unsigned int r,
                                                              const unsigned int s ) const
  {
    return _Mr_Ms_to_the_one_fourth[r][s];
  }
    
  
  template<class CoeffType>
  inline
  CoeffType WilkeMixture<CoeffType>::denominator( const unsigned int r,
                                                  const unsigned int s ) const
  {
    return _denom[r][s];
  }

  template<class CoeffType>
  inline
  const ChemicalMixture<CoeffType>& WilkeMixture<CoeffType>::chem_mixture() const
  {
    return _chem_mixture;
  }

} // end namespace Antioch

#endif // ANTIOCH_WILKE_MIXTURE_H
