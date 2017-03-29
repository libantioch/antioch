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

#ifndef ANTIOCH_MICRO_THERMO_BASE_H
#define ANTIOCH_MICRO_THERMO_BASE_H

// Antioch
#include "antioch/chemical_mixture.h"

namespace Antioch
{

  template<typename CoeffType>
  class MicroThermoBase
  {
  public:

    MicroThermoBase( const ChemicalMixture<CoeffType>& chem_mixture )
      : _chem_mixture(chem_mixture)
    {}

    //! Pure virtual destructor
    /*! This is pure virtual to force this object to be abstract. */
    virtual ~MicroThermoBase() =0;

  protected:

    const ChemicalMixture<CoeffType> & _chem_mixture;

  private:

    //! Default constructor
    /*! Private to force to user to supply a ChemicalMixture object.*/
    MicroThermoBase();
  };

  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType>
  inline
  MicroThermoBase<CoeffType>::~MicroThermoBase(){}

} // end namespace Antioch

#endif // ANTIOCH_MICRO_THERMO_BASE_H
