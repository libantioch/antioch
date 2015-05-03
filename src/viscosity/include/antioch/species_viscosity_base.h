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

#ifndef ANTIOCH_SPECIES_VISCOSITY_BASE_H
#define ANTIOCH_SPECIES_VISCOSITY_BASE_H

// C++
#include <vector>
#include <ostream>

namespace Antioch
{
  //! Base class for species viscosity models
  /*! We use the curiously recurring template pattern to enforce
      the interface that subclasses must adhere in order to
      ulimately be used in the MixtureViscosity class. Subclasses
      must implement:
         -# op_impl --- this should implement operator()
         -# reset_coeffs_impl --- should implement reset_coeffs
         -# print_impl --- should implement print
  */
  template<typename Subclass, typename CoeffType>
  class SpeciesViscosityBase
  {
  public:

    SpeciesViscosityBase(){};

    virtual ~SpeciesViscosityBase(){};

    //! Evaluates viscosity at temperature T
    template <typename StateType>
    StateType operator()( const StateType& T ) const;

    //! Resets coefficients associated with the viscosity model
    void reset_coeffs( const std::vector<CoeffType>& coeffs );

    //! Formatted print, by default to \p std::cout
    void print(std::ostream& os = std::cout) const;

    //! Formatted print.
    friend std::ostream& operator<<(std::ostream& os, const Subclass& mu)
    {
      mu.print(os);
      return os;
    }
  };

  template<typename Subclass, typename CoeffType>
  template <typename StateType>
  inline
  StateType SpeciesViscosityBase<Subclass,CoeffType>::operator()( const StateType& T ) const
  {
    return static_cast<const Subclass*>(this)->op_impl(T);
  }

  template<typename Subclass, typename CoeffType>
  void SpeciesViscosityBase<Subclass,CoeffType>::reset_coeffs( const std::vector<CoeffType>& coeffs )
  {
    static_cast<Subclass*>(this)->reset_coeffs_impl(coeffs);
  }

  template<typename Subclass, typename CoeffType>
  void SpeciesViscosityBase<Subclass,CoeffType>::print(std::ostream& os) const
  {
    static_cast<const Subclass*>(this)->print_impl(os);
  }

} // end namespace Antioch

#endif // ANTIOCH_SPECIES_VISCOSITY_BASE_H
