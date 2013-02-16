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

#include "antioch/chemical_species.h"

namespace Antioch
{
  template<class NumericType>
  ChemicalSpecies<NumericType>::ChemicalSpecies( const std::string &name,
						 const NumericType mol_wght,
						 const NumericType h_form,
						 const NumericType n_tr_dofs,
						 const int charge )
    : _name      (name),
      _mol_wght  (mol_wght),
      _R         (Constants::R_universal/mol_wght),
      _h_form    (h_form),
      _n_tr_dofs (n_tr_dofs),
      _charge    (charge)
  {
    return;
  }

  template<class NumericType>
  ChemicalSpecies<NumericType>::ChemicalSpecies()
    : _name("Err"), 
      _mol_wght(0), 
      _R(0),
      _h_form(0), 
      _n_tr_dofs(0), 
      _charge(0)
  {
    antioch_error();
    return;
  }

  template<class NumericType>
  ChemicalSpecies<NumericType>::~ChemicalSpecies()
  {
    return;
  }

  template<class NumericType>
  void ChemicalSpecies<NumericType>::print (std::ostream &os) const
  {
    os << " -----------------------------\n"
       << "| Species " << this->species() << '\n'
       << " -----------------------------\n"
       << std::scientific
       << "  Mol Wgt = " << this->molar_mass() << '\n'
       << "  R       = " << this->gas_constant()     << '\n'
       << "  h0      = " << this->formation_enthalpy() << '\n'
       << "  n_tr    = " << this->n_tr_dofs() << '\n'
       << "  charge  = " << this->charge() << '\n';

    for (unsigned int l=0; l<this->theta_v().size(); l++)
      os << "  theta_v_" << l << " = " << this->theta_v()[l]
	 << ", ndg = " << this->ndg_v()[l] << "\n";
    
    for (unsigned int l=0; l<this->theta_e().size(); l++)
      os << "  theta_e_" << l << " = " << this->theta_e()[l]
	 << ", ndg = " << this->ndg_e()[l] << "\n";
    
    os << '\n';
    
    return;
  }

  /* ------------------------- Instantiate ------------------------- */
  template class ChemicalSpecies<double>;

} // end namespace Antioch
