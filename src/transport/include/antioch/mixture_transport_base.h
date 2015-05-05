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

#ifndef ANTIOCH_MIXTURE_TRANSPORT_BASE_H
#define ANTIOCH_MIXTURE_TRANSPORT_BASE_H

#include "antioch/transport_mixture.h"

namespace Antioch
{
  template<typename ThermoEvaluator, class CoeffType=double>
  class MixtureTransportBase
  {
  public:

    MixtureTransportBase( const TransportMixture<ThermoEvaluator,CoeffType>& transport_mixture );

    virtual ~MixtureTransportBase(){};

    const ChemicalMixture<CoeffType>& chemical_mixture() const;

    const TransportMixture<ThermoEvaluator,CoeffType>& transport_mixture() const;

    const TransportMixture<ThermoEvaluator,CoeffType>& mixture() const;

  protected:

    const TransportMixture<ThermoEvaluator,CoeffType>& _transport_mixture;

  private:

    MixtureTransportBase();

  };

  template<typename ThermoEvaluator, class CoeffType>
  MixtureTransportBase<ThermoEvaluator,CoeffType>::MixtureTransportBase( const TransportMixture<ThermoEvaluator,CoeffType>& transport_mixture )
    : _transport_mixture(transport_mixture)
  {}

  template<typename ThermoEvaluator, class CoeffType>
  inline
  const ChemicalMixture<CoeffType>& MixtureTransportBase<ThermoEvaluator,CoeffType>::chemical_mixture() const
  {
    return _transport_mixture.chemical_mixture();
  }

  template<typename ThermoEvaluator, class CoeffType>
  inline
  const TransportMixture<ThermoEvaluator,CoeffType>& MixtureTransportBase<ThermoEvaluator,CoeffType>::transport_mixture() const
  {
    return _transport_mixture;
  }

  template<typename ThermoEvaluator, class CoeffType>
  inline
  const TransportMixture<ThermoEvaluator,CoeffType>& MixtureTransportBase<ThermoEvaluator,CoeffType>::mixture() const
  {
    return _transport_mixture;
  }

} // end namespace Antioch

#endif // ANTIOCH_MIXTURE_TRANSPORT_BASE_H
