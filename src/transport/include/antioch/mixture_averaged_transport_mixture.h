//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014-2016 Paul T. Bauman, Benjamin S. Kirk,
//                         Sylvain Plessis, Roy H. Stonger
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

#ifndef ANTIOCH_WILKE_TRANSPORT_MIXTURE_H
#define ANTIOCH_WILKE_TRANSPORT_MIXTURE_H

// C++
#include <vector>

namespace Antioch
{
  // Forward declarations
  template <typename CoeffType>
  class ChemicalMixture;

  template <typename CoeffType>
  class TransportMixture;

  //! Mixture object for MixtureAveragedTransport model
  /*! This object is meant to live for the life of the program and contains data that is reused
   *  by the MixtureAveragedTransportEvaluator object. In particular, there are two terms for
   *  species indices r and s:
   *
   * \f[  \left(\frac{M_r}{M_s}\right)^{1/4} \f]
   *
   * \f[ \sqrt{8\left( 1 + \frac{M_r}{M_s} \right)} \f]
   *
   * These terms appear in the mixing formulae used in MixtureAveragedTransportEvaluator.
   */
  template<class CoeffType = double>
  class MixtureAveragedTransportMixture
  {
  public:

    MixtureAveragedTransportMixture( const TransportMixture<CoeffType> & mixture);

    ~MixtureAveragedTransportMixture(){};

    //! \f[ \left(\frac{M_r}{M_s}\right)^{1/4} \f]
    CoeffType Mr_Ms_to_the_one_fourth( const unsigned int r,
                                       const unsigned int s ) const;

    //! \f[ \sqrt{8\left( 1 + \frac{M_r}{M_s} \right)} \f]
    CoeffType denominator( const unsigned int r,
                           const unsigned int s ) const;

    //! chemical mixture, mostly for backward compatibility
    const ChemicalMixture<CoeffType>& chem_mixture() const;

    //! transport mixture
    const TransportMixture<CoeffType> & transport_mixture() const;

  protected:

    const TransportMixture<CoeffType> & _mixture;

    //! Cache for numerator term
    /*! \todo We should use a more efficient data structure */
    std::vector<std::vector<CoeffType> > _Mr_Ms_to_the_one_fourth;

    //! Cache for denominator term
    /*! \todo We should use a more efficient data structure */
    std::vector<std::vector<CoeffType> > _denom;

  };

  template<class CoeffType>
  MixtureAveragedTransportMixture<CoeffType>::MixtureAveragedTransportMixture( const TransportMixture<CoeffType>& mixture)
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
  }

  template<class CoeffType>
  inline
  CoeffType MixtureAveragedTransportMixture<CoeffType>::Mr_Ms_to_the_one_fourth( const unsigned int r,
                                                                                 const unsigned int s ) const
  {
    return _Mr_Ms_to_the_one_fourth[r][s];
  }


  template<class CoeffType>
  inline
  CoeffType MixtureAveragedTransportMixture<CoeffType>::denominator( const unsigned int r,
                                                                     const unsigned int s ) const
  {
    return _denom[r][s];
  }

  template<class CoeffType>
  inline
  const ChemicalMixture<CoeffType>& MixtureAveragedTransportMixture<CoeffType>::chem_mixture() const
  {
    return _mixture.chemical_mixture();
  }

  template<class CoeffType>
  inline
  const TransportMixture<CoeffType> & MixtureAveragedTransportMixture<CoeffType>::transport_mixture() const
  {
    return _mixture;
  }

} // end namespace Antioch

#endif // ANTIOCH_WILKE_TRANSPORT_MIXTURE_H
