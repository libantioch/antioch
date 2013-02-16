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

#ifndef ANTIOCH_CEA_THERMO_H
#define ANTIOCH_CEA_THERMO_H

// C++
#include <iomanip>
#include <vector>

namespace Antioch
{
  // Forward declarations
  template<class NumericType>
  class ChemicalMixture;

  template<class NumericType>
  class CEACurveFit;

  template<class NumericType>
  class CEAThermodynamics
  {
  public:

    CEAThermodynamics( const ChemicalMixture<NumericType>& chem_mixture );

    //! Destructor
    /*! virtual so this can be subclassed by the user. */
    virtual ~CEAThermodynamics();

    NumericType cp( NumericType T, unsigned int species ) const;

    NumericType cp( NumericType T, const std::vector<NumericType>& mass_fractions ) const;

    NumericType cv( NumericType T, unsigned int species ) const;

    NumericType cv( NumericType T, const std::vector<NumericType>& mass_fractions ) const;

    NumericType h( NumericType T, unsigned int species ) const;

    void h( NumericType T, std::vector<NumericType>& h ) const;

    NumericType h_RT_minus_s_R( NumericType T, unsigned int species ) const;

    void h_RT_minus_s_R( NumericType T, std::vector<NumericType>& h_RT_minus_s_R ) const;

    NumericType cp_over_R( NumericType T, unsigned int species ) const;

    NumericType h_over_RT( NumericType T, unsigned int species ) const;

    NumericType s_over_R( NumericType T, unsigned int species ) const;

  protected:

    void read_thermodynamic_table();

    void read_thermodynamic_table( std::istream& in );

    const ChemicalMixture<NumericType>& _chem_mixture;

    std::vector<CEACurveFit<NumericType>* > _species_curve_fits;

    std::vector<NumericType> _cp_at_200p1;

  private:
    
    //! Default constructor
    /*! Private to force to user to supply a ChemicalMixture object.*/
    CEAThermodynamics();

  };
  
  /* ------------------------- Inline Functions -------------------------*/
  
  template<class NumericType>
  inline
  NumericType CEAThermodynamics<NumericType>::cv( NumericType T, unsigned int species ) const
  {
    return this->cp(T,species) - _chem_mixture.R(species);
  }

  template<class NumericType>
  inline
  NumericType CEAThermodynamics<NumericType>::cv( NumericType T,
						  const std::vector<NumericType>& mass_fractions ) const
  {
    return this->cp(T,mass_fractions) - _chem_mixture.R(mass_fractions);
  }

} // end namespace Antioch

#endif //ANTIOCH_CEA_THERMO_H
