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

#ifndef ANTIOCH_BLOTTNER_PARSING_H
#define ANTIOCH_BLOTTNER_PARSING_H

// Antioch
#include "antioch/blottner_viscosity.h"
#include "antioch/mixture_viscosity.h"

// C++
#include <iostream>

namespace Antioch
{
  template<class NumericType>
  void read_blottner_data_ascii( MixtureViscosity<BlottnerViscosity<NumericType>,NumericType >& mu,
				 std::istream &in );

  template<class NumericType>
  void read_blottner_data_ascii_default( MixtureViscosity<BlottnerViscosity<NumericType>,NumericType>& mu );


  /* ------------------------- Inline Functions -------------------------*/
  template<class NumericType>
  inline
  void read_blottner_data_ascii( MixtureViscosity<BlottnerViscosity<NumericType>,NumericType >& mu,
				 std::istream &in )
  {
    // skip the header
    skip_comment_lines(in, '#');

    std::string name;
    NumericType a, b, c;

    while (in.good())
      {
        in >> name; // Species Name
        in >> a;    // 
        in >> b;    // 
        in >> c;    // 
        
        // If we are still good, we have a valid set of transport
        // data for this species. Otherwise, we read past end-of-file 
        // in the section above
        if (in.good())
          {
	    const ChemicalMixture<NumericType>& chem_mixture = mu.chemical_mixture();
	    
	    // Check if this is a species we want.
	    if( chem_mixture.active_species_name_map().find(name) !=
		chem_mixture.active_species_name_map().end() )
	      {
		// Pack up coefficients
		std::vector<NumericType> coeffs(3);
		coeffs[0] = a;
		coeffs[1] = b;
		coeffs[2] = c;
		mu.add(name, coeffs);
	      }
          }

      }

    return;
  }
  

  template<class NumericType>
  inline
  void read_blottner_data_ascii_default( MixtureViscosity<BlottnerViscosity<NumericType>,NumericType>& mu )
  {
    static const std::string
      default_blottner_transport_data
      ("#-----------------------------------------------------------------------------\n"
       "# Coefficients for Blottner viscosity model\n"
       "# Form of the fit:\n"
       "# \n"
       "# mu = 0.1*exp[A*log(T)*log(T) + B*log(T) + C]\n"
       "# \n"
       "# References:\n"
       "# \n"
       "# Air  -- Identical to N2 fit\n"
       "# N    -- Sandia Report SC-RR-70-754\n"
       "# N2   -- Sandia Report SC-RR-70-754\n"
       "# CPN2 -- Identical to N2 fit\n"
       "# NO   -- Sandia Report SC-RR-70-754\n"
       "# O    -- Sandia Report SC-RR-70-754\n"
       "# O2   -- Sandia Report SC-RR-70-754\n"
       "# C    -- AIAA-1997-2474\n"
       "# C2   -- AIAA-1997-2474\n"
       "# C3   -- AIAA-1997-2474\n"
       "# C2H  -- wild-ass guess: identical to HCN fit\n"
       "# CN   -- AIAA-1997-2474\n"
       "# CO   -- AIAA-1997-2474\n"
       "# CO2  -- AIAA-1997-2474\n"
       "# HCN  -- AIAA-1997-2474\n"
       "# H    -- AIAA-1997-2474\n"
       "# H2   -- AIAA-1997-2474\n"
       "# e    -- Sandia Report SC-RR-70-754\n"
       "\n"
       "Air     2.68142000000e-02  3.17783800000e-01 -1.13155513000e+01\n"
       "CPAir   2.68142000000e-02  3.17783800000e-01 -1.13155513000e+01\n"
       "N       1.15572000000e-02  6.03167900000e-01 -1.24327495000e+01\n"
       "N2      2.68142000000e-02  3.17783800000e-01 -1.13155513000e+01\n"
       "CPN2    2.68142000000e-02  3.17783800000e-01 -1.13155513000e+01\n"
       "NO      4.36378000000e-02 -3.35511000000e-02 -9.57674300000e+00\n"
       "O       2.03144000000e-02  4.29440400000e-01 -1.16031403000e+01\n"
       "O2      4.49290000000e-02 -8.26158000000e-02 -9.20194750000e+00\n"
       "C       -8.3285e-3         0.7703240         -12.7378000\n"
       "C2      -8.4311e-3         0.7876060         -13.0268000\n"
       "C3      -8.4312e-3         0.7876090         -12.8240000\n"
       "C2H     -2.4241e-2         1.0946550         -14.5835500\n"
       "CN      -8.3811e-3         0.7860330         -12.9406000\n"
       "CO      -0.019527394       1.013295          -13.97873\n"
       "CO2     -0.019527387       1.047818          -14.32212\n"
       "HCN     -2.4241e-2         1.0946550         -14.5835500\n"
       "H       -8.3912e-3         0.7743270         -13.6653000\n"
       "H2      -8.3346e-3         0.7815380         -13.5351000\n"
       "e       0.00000000000e+00  0.00000000000e+00 -1.16031403000e+01\n");

    std::istringstream buf(default_blottner_transport_data);

    read_blottner_data_ascii( mu, buf );

    return;
  }


} // end namespace Antioch

#endif // ANTIOCH_BLOTTNER_PARSING_H
