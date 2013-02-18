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
#include <cmath>

// Antioch
#include "antioch/chemical_mixture.h"

namespace Antioch
{
  // Forward declarations
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

    class Cache 
    {
    public:
      const NumericType &T;
      NumericType T2;
      NumericType T3;
      NumericType T4;
      NumericType lnT;
      
      Cache(const NumericType &T_in) 
	: T(T_in)
      {
	T2  = T*T;
	T3  = T2*T;
	T4  = T2*T2;
	lnT = std::log(T);	
      }
    private:
      Cache();      
    };
    
    void add_curve_fit( const std::string& species_name, const std::vector<NumericType>& coeffs );

    //! Checks that curve fits have been specified for all species in the mixture.
    bool check() const;

    NumericType cp( const Cache &cache, unsigned int species ) const;

    NumericType cp( const Cache &cache, const std::vector<NumericType>& mass_fractions ) const;

    NumericType cv( const Cache &cache, unsigned int species ) const;

    NumericType cv( const Cache &cache, const std::vector<NumericType>& mass_fractions ) const;

    NumericType h( const Cache &cache, unsigned int species ) const;

    void h( const Cache &cache, std::vector<NumericType>& h ) const;

    NumericType h_RT_minus_s_R( const Cache &cache, unsigned int species ) const;

    void h_RT_minus_s_R( const Cache &cache, std::vector<NumericType>& h_RT_minus_s_R ) const;

    NumericType cp_over_R( const Cache &cache, unsigned int species ) const;

    NumericType h_over_RT( const Cache &cache, unsigned int species ) const;

    NumericType s_over_R( const Cache &cache, unsigned int species ) const;

    const ChemicalMixture<NumericType>& chemical_mixture() const;

  protected:

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
  NumericType CEAThermodynamics<NumericType>::cv( const Cache &cache, unsigned int species ) const
  {
    return this->cp(cache,species) - _chem_mixture.R(species);
  }

  template<class NumericType>
  inline
  NumericType CEAThermodynamics<NumericType>::cv( const Cache &cache,
						  const std::vector<NumericType>& mass_fractions ) const
  {
    return this->cp(cache,mass_fractions) - _chem_mixture.R(mass_fractions);
  }

  template<class NumericType>
  inline
  const ChemicalMixture<NumericType>& CEAThermodynamics<NumericType>::chemical_mixture() const
  {
    return _chem_mixture;
  }

} // end namespace Antioch

#endif //ANTIOCH_CEA_THERMO_H
