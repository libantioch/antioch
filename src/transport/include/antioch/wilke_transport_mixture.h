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

#ifndef ANTIOCH_WILKE_TRANSPORT_MIXTURE_H
#define ANTIOCH_WILKE_TRANSPORT_MIXTURE_H

// C++
#include <vector>

namespace Antioch
{

  // forward declaration to keep the call to 
  // the ChemicalMixture object
  template <typename CoeffType>
  class ChemicalMixture;

  template <typename ThermoEval, typename CoeffType>
  class TransportMixture;

  template<class ThermoEvaluator,         // ThermoHandler
           class CoeffType = double>      // type
  class WilkeTransportMixture
  {
  public:

    WilkeTransportMixture( const TransportMixture<ThermoEvaluator,CoeffType> & mixture);
    ~WilkeTransportMixture();

    CoeffType Mr_Ms_to_the_one_fourth( const unsigned int r,
                                       const unsigned int s ) const;
    
    CoeffType denominator( const unsigned int r,
                           const unsigned int s ) const;

    //! chemical mixture, mostly for backward compatibility
    const ChemicalMixture<CoeffType>& chem_mixture() const;

    //! transport mixture
    const TransportMixture<ThermoEvaluator,CoeffType> & transport_mixture() const;

    //! const ref to thermo
    const ThermoEvaluator & thermo_evaluator() const;

  protected:

    const TransportMixture<ThermoEvaluator,CoeffType> & _mixture;

    //! Cache for numerator term
    /*! \todo We should use a more efficient data structure */
    std::vector<std::vector<CoeffType> > _Mr_Ms_to_the_one_fourth;
    
    //! Cache for denominator term
    /*! \todo We should use a more efficient data structure */
    std::vector<std::vector<CoeffType> > _denom;

  };

  template<class ThermoEvaluator, class CoeffType>
  WilkeTransportMixture<ThermoEvaluator,CoeffType>::WilkeTransportMixture( const TransportMixture<ThermoEvaluator,CoeffType>& mixture)
    : _mixture(mixture),
      _Mr_Ms_to_the_one_fourth(mixture.n_species()),
      _denom(mixture.n_species())
  {
    using std::pow;

    for( unsigned int r = 0; r < mixture.n_species(); r++ )
      {
        _Mr_Ms_to_the_one_fourth[r].resize(mixture.n_species());
        _denom[r].resize(mixture.n_species());

        for( unsigned int s = 0; s < mixture.n_species(); s++ )
          {
            const CoeffType Mr = mixture.chemical_mixture().M(r);
            const CoeffType Ms = mixture.chemical_mixture().M(s);

            _Mr_Ms_to_the_one_fourth[r][s] = pow( Mr/Ms, CoeffType(0.25) );
            _denom[r][s] = std::sqrt(8.0*(1.0+Ms/Mr));
          }
      }

    return;
  }

  template<class ThermoEvaluator, class CoeffType>
  WilkeTransportMixture<ThermoEvaluator,CoeffType>::~WilkeTransportMixture()
  {
    return;
  }

  
  template<class ThermoEvaluator, class CoeffType>
  inline
  CoeffType WilkeTransportMixture<ThermoEvaluator,CoeffType>::Mr_Ms_to_the_one_fourth( const unsigned int r,
                                                              const unsigned int s ) const
  {
    return _Mr_Ms_to_the_one_fourth[r][s];
  }
    
  
  template<class ThermoEvaluator, class CoeffType>
  inline
  CoeffType WilkeTransportMixture<ThermoEvaluator,CoeffType>::denominator( const unsigned int r,
                                                  const unsigned int s ) const
  {
    return _denom[r][s];
  }

  template<class ThermoEvaluator, class CoeffType>
  inline
  const ChemicalMixture<CoeffType>& WilkeTransportMixture<ThermoEvaluator,CoeffType>::chem_mixture() const
  {
    return _mixture.chemical_mixture();
  }

  template<class ThermoEvaluator, class CoeffType>
  inline
  const TransportMixture<ThermoEvaluator,CoeffType> & WilkeTransportMixture<ThermoEvaluator,CoeffType>::transport_mixture() const
  {
    return _mixture;
  }

  template<class ThermoEvaluator, class CoeffType>
  inline
  const ThermoEvaluator & WilkeTransportMixture<ThermoEvaluator,CoeffType>::thermo_evaluator() const
  {
      return _mixture.thermo();
  }

} // end namespace Antioch

#endif // ANTIOCH_WILKE_TRANSPORT_MIXTURE_H
