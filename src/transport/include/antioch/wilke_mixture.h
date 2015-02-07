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

#ifndef ANTIOCH_WILKE_MIXTURE_H
#define ANTIOCH_WILKE_MIXTURE_H

// C++
#include <vector>

namespace Antioch
{

  // forward declaration to keep the call to 
  // the ChemicalMixture object
  template <typename CoeffType>
  class ChemicalMixture;

  /*
      This somewhat ridicule thermodynamics
      construct is because StatMechThermodynamics
      can't give cp, and NASA9 can't give
      the vibrational, translational, etc.,
      parts of this same quantity.
      Notice that the internal R used is the
      massic one, thus Cp is in J/kg/K

      \todo make Wilke rules be computed here?
  */
  template<class Mixture, class ThermoEvaluator,                      // mixture + thermo
           class CoeffType = double>                                  // type
  class WilkeMixture
  {
  public:

    WilkeMixture( const Mixture & mixture,
                  const ThermoEvaluator & thermo_eval );
    ~WilkeMixture();

    CoeffType Mr_Ms_to_the_one_fourth( const unsigned int r,
                                       const unsigned int s ) const;
    
    CoeffType denominator( const unsigned int r,
                           const unsigned int s ) const;

    //! chemical mixture, mostly for backward compatibility
    const ChemicalMixture<CoeffType>& chem_mixture() const;

    //! transport mixture
    const Mixture & transport_mixture() const;  // contains the macro thermo for species

    //! const ref to thermo evaluator, internal computations only (vib, rot, trans)
    const ThermoEvaluator & thermo_evaluator() const; // contains the micro thermo

  protected:

    const Mixture         & _mixture;

    const ThermoEvaluator & _thermo_evaluator;

    //! Cache for numerator term
    /*! \todo We should use a more efficient data structure */
    std::vector<std::vector<CoeffType> > _Mr_Ms_to_the_one_fourth;
    
    //! Cache for denominator term
    /*! \todo We should use a more efficient data structure */
    std::vector<std::vector<CoeffType> > _denom;

  };

  template<class Mixture, class ThermoEvaluator, class CoeffType>
  WilkeMixture<Mixture,ThermoEvaluator,CoeffType>::WilkeMixture( const Mixture& mixture, const ThermoEvaluator & thermo_eval )
    : _mixture(mixture),
      _thermo_evaluator(thermo_eval),
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

  template<class Mixture, class ThermoEvaluator, class CoeffType>
  WilkeMixture<Mixture,ThermoEvaluator,CoeffType>::~WilkeMixture()
  {
    return;
  }

  
  template<class Mixture, class ThermoEvaluator, class CoeffType>
  inline
  CoeffType WilkeMixture<Mixture,ThermoEvaluator,CoeffType>::Mr_Ms_to_the_one_fourth( const unsigned int r,
                                                              const unsigned int s ) const
  {
    return _Mr_Ms_to_the_one_fourth[r][s];
  }
    
  
  template<class Mixture, class ThermoEvaluator, class CoeffType>
  inline
  CoeffType WilkeMixture<Mixture,ThermoEvaluator,CoeffType>::denominator( const unsigned int r,
                                                  const unsigned int s ) const
  {
    return _denom[r][s];
  }

  template<class Mixture, class ThermoEvaluator, class CoeffType>
  inline
  const ChemicalMixture<CoeffType>& WilkeMixture<Mixture,ThermoEvaluator,CoeffType>::chem_mixture() const
  {
    return _mixture.chemical_mixture();
  }

  template<class Mixture, class ThermoEvaluator, class CoeffType>
  inline
  const Mixture & WilkeMixture<Mixture,ThermoEvaluator,CoeffType>::transport_mixture() const
  {
    return _mixture;
  }

  template<class Mixture, class ThermoEvaluator, class CoeffType>
  inline
  const ThermoEvaluator & WilkeMixture<Mixture,ThermoEvaluator,CoeffType>::thermo_evaluator() const
  {
      return _thermo_evaluator;
  }

} // end namespace Antioch

#endif // ANTIOCH_WILKE_MIXTURE_H
