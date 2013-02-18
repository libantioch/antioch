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

// C++
#include <sstream>

// These functions
#include "antioch/species_ascii_parsing.h"

// Antioch
#include "antioch/input_utils.h"
#include "antioch/chemical_mixture.h"

namespace Antioch
{

  template<class NumericType>
  void read_species_data_ascii( ChemicalMixture<NumericType>& chem_mixture,
				std::istream& in )
  {
    skip_comment_lines(in, '#');

    std::string name;
    NumericType mol_wght, h_form, n_tr_dofs;
    int charge;

    while (in.good())
      {
	in >> name;      // Species Name
	in >> mol_wght;  // molecular weight (kg/kmol)
	in >> h_form;    // heat of formation at Ok (J/kg)
	in >> n_tr_dofs; // number of translational/rotational DOFs
	in >> charge;    // charge number
	
	// If we are still good, we have a valid set of thermodynamic
	// data for this species. Otherwise, we read past end-of-file 
	// in the section above
	if (in.good())
	  {
	    // If we do not have an enum for this species, that is going to cause
	    // problems down the line.
	    if (!chem_mixture.species_name_map().count(name))
	      {
		std::cerr << "ERROR: Unexpected species " << name 
			  << " encountered while parsing thermodynamic table!"
			  << std::endl;
		antioch_error();
	      }

	    /* Only add a ChemicalSpecies if the user supplied it to
	       the ChemicalMixture. */
	    Species species = chem_mixture.species_name_map().find(name)->second;

	    // using default comparison:
	    std::vector<Species>::const_iterator it = std::search_n( chem_mixture.species_list().begin(), 
								     chem_mixture.species_list().end(),
								     1, species);
	    if( it != chem_mixture.species_list().end() )
	      {
		unsigned int index = static_cast<unsigned int>(it - chem_mixture.species_list().begin());
		chem_mixture.add_species( index, name, mol_wght, h_form, n_tr_dofs, charge );
	      }

	  }
      }

    return;
  }

  template<class NumericType>
  void read_species_data_ascii_default( ChemicalMixture<NumericType>& chem_mixture )
  {
    /*!
     * Default species data.  Includes molecular weights, 
     * heats of formation (J/kg), nominal translational/rotational
     * degrees of freedom (cv=c_tr*R*T where c_tr=5/2 for diatomics
     * and c_tr=3/2 for atoms, and charge number.)
     * Originally taken from FIN-S.
     */
    static const std::string
      default_species_chem_data
      ("#===========================================================================\n"
       "# LEGEND\n"
       "#===========================================================================\n"
       "#\n"    
       "#    Species 	-- Species name\n"
       "#    Mol. Wt.	-- Molecular weight (kg/kmol)\n"
       "#    Hform	-- Formation enthalpy (J/kg @ 0K)\n"
       "#    cfs	-- Nominal T-R Degrees of freedom (cv = cfs*k*T)\n"
       "#    zns	-- Charge number\n"
       "#\n"
       "#[PB]: I'm unsure about the cfs for HO2 and H2O2."
       "#    Spec.     Mol. Wt.     Hform (0K),  cfs,   zns\n"
       "\n"    
       "Air	28.96000	  0.000000000000 	2.5 	 0\n"
       "CPAir	28.96000	  0.000000000000 	2.5 	 0\n"
       "Ar	39.94400	  0.000000000000 	1.5 	 0\n"
       "Ar+	39.94345	  3.8068120000e7 	1.5 	 1\n"
       "C	12.01100	  5.9211889000e7 	1.5 	 0\n"
       "C+	12.01045	  1.4967366000e8 	1.5 	 1\n"
       "C2	24.02200	  3.4234785000e7 	2.5 	 0\n"
       "C2H	25.03000	  2.2572910000e7 	2.5 	 0\n"
       "C2H2	26.03800	  8.7548580000e6 	2.5 	 0\n"
       "C3      36.03300	  2.3062193000e7 	2.5 	 0\n"
       "CF	31.00940	  9.9617550000e6 	2.5 	 0\n"
       "CF2	50.00780	 -2.5187600000e6 	3.0 	 0\n"
       "CF3	69.00620	 -7.1992350000e6 	3.0 	 0\n"
       "CF4	88.00460	 -1.0258770000e7 	3.0 	 0\n"
       "CH	13.01900	  4.5627850000e7 	2.5 	 0\n"
       "CH2	14.02700	  2.7803520000e7 	3.0 	 0\n"
       "CH3	15.03500	  9.9559030000e6 	3.0 	 0\n"
       "CH4	16.04300	 -4.1532130000e6 	3.0 	 0\n"
       "Cl      35.45300	  3.3740400000e6 	1.5 	 0\n"
       "Cl2     70.90600	  0.000000000000 	2.5 	 0\n"
       "CN	26.01900	  1.6795420000e7 	2.5 	 0\n"
       "CN+	26.01845	  6.8835800000e7 	2.5 	 1\n"
       "CO	28.01100	 -4.0630824000e6 	2.5 	 0\n"
       "CO+	28.01045	  4.4200904000e7 	2.5 	 1\n"
       "CO2	44.01100	 -8.9328800000e6 	2.5 	 0\n"
       "F	18.99840	  4.0423050000e6 	1.5 	 0\n"
       "F2	37.99680	  0.000000000000 	2.5 	 0\n"
       "H	 1.00800	  2.1432040000e8 	1.5 	 0\n"
       "H+	 1.00745	  1.5167840000e9 	1.5 	 1\n"
       "H2	 2.01600	  0.000000000000 	2.5 	 0\n"
       "H2+	 2.01545	  7.3847530000e8 	2.5 	 1\n"
       "H2O	18.01600	 -1.3261710000e7 	3.0 	 0\n"
       "H2O2    34.01475         -3.8186669018e6        3.0      0\n"
       "HCl     36.46100	 -2.5266800000e6 	2.5 	 0\n"
       "HCN	27.02700	  4.8982130000e6 	2.5 	 0\n"
       "He	 4.00300	  0.000000000000 	1.5 	 0\n"
       "He+	 4.00245	  5.9271800000e8 	1.5 	 1\n"
       "HO2     33.00681          4.5239149133e5        3.0      0\n"
       "N	14.00800	  3.3621610000e7 	1.5 	 0\n"
       "Ne	20.17900	  0.000000000000 	1.5 	 0\n"
       "N+	14.00745	  1.3400000000e8 	1.5 	 1\n"
       "N2	28.01600	  0.000000000000 	2.5 	 0\n"
       "CPN2	28.01600	  0.000000000000 	2.5 	 0\n"
       "N2+	28.01545	  5.3700000000e7 	2.5 	 1\n"
       "NCO	42.01900	  4.2124000000e6 	2.5 	 0\n"
       "NH	15.01600	  2.3867900000e7 	2.5 	 0\n"
       "NH+	15.01545	  1.1050000000e8 	2.5 	 1\n"
       "NH2	16.02400	  1.2036000000e7 	3.0 	 0\n"
       "NH3	17.03200	 -2.2866370000e6 	3.0 	 0\n"
       "NO	30.00800	  2.9961230000e6 	2.5 	 0\n"
       "NO+	30.00745	  3.2834800000e7 	2.5 	 1\n"
       "NO2	46.00800	  8.0420800000e5 	3.0 	 0\n"
       "O	16.00000	  1.5420000000e7 	1.5 	 0\n"
       "O+	15.99945	  9.7560000000e7 	1.5 	 1\n"
       "O2	32.00000	  0.000000000000 	2.5 	 0\n"
       "O2+	31.99945	  3.6370000000e7 	2.5 	 1\n"
       "OH	17.00800	  2.2995060000e6 	2.5 	 0\n"
       "Si	28.08550	  1.5868220000e7 	1.5 	 0\n"
       "SiO	44.08550	 -2.2683200000e6 	2.5 	 0\n"
       "e	 0.00055	  0.000000000000 	1.5 	-1\n");

    std::istringstream buf(default_species_chem_data);
    read_species_data_ascii(chem_mixture,buf);

    return;
  }

  /* ------------------------- Instantiate ------------------------- */
  template void read_species_data_ascii<double>( ChemicalMixture<double>&, std::istream & );
  template void read_species_data_ascii_default<double>( ChemicalMixture<double>& );

} // namespace Antioch
