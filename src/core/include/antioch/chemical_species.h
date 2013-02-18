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

#ifndef ANTIOCH_CHEMICAL_SPECIES_H
#define ANTIOCH_CHEMICAL_SPECIES_H

// C++ includes
#include <string>
#include <vector>
#include <iostream>

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/physical_constants.h"

namespace Antioch
{
  //! Class to encapsulate data for each chemical species
  /*!
   * This class is designed to store information relevant to a chemical species.
   * All the data stored is constant for each species, so we store const for each
   * variable. The idea is that this will be placed inside ChemicalMixture, which
   * will be a singleton. This is stolen from the FIN-S class SpeciesChemistry.
   */
  template<class NumericType>
  class ChemicalSpecies
  {
  public:
    
    //! Constrctor
    ChemicalSpecies( const std::string &name, 
		     const NumericType mol_wght,
		     const NumericType h_form,
		     const NumericType n_tr_dofs,
		     const int         charge );

    //! Default constructor.
    /*!
     * This is technically required for any
     * std::map value type (or operator[] breaks, at least).  But,
     * we never actually want to create a SpeciesChemistry
     * implicitly, so we throw an error if this is ever used.
     */      
    ChemicalSpecies();

    //! Destructor
    ~ChemicalSpecies();

    //! Returns a descriptive name for this species.
    const std::string& species() const;

    //!Returns the molar mass in (g/mol) or (kg/kmol).
    NumericType molar_mass() const;

    //! Returns the species ideal gas constant in [J/kg-K]
    /*!
     * \f$ R \equiv \frac{\hat{R}}{M} \f$ where
     * \f$ R\f$ is the universal gas constant and
     * \f$ M \f$ is the species molar mass.
     */
    NumericType gas_constant() const;

    //! Returns formation enthalpy in units of [J/kg]
    NumericType formation_enthalpy() const;
    
    //! Returns number of translational degrees of freedom
    NumericType n_tr_dofs() const;
    
    //! Returns electrical charge number
    int charge() const;

    //! Returns true if the chemical species has vibrational degrees 
    bool has_vibrational_modes() const;

    //! Returns true if the chemical species has vibrational degrees of freedom.
    unsigned int n_vibrational_modes() const;

    //!Characteristic vibrational temperature [K].
    const std::vector<NumericType>& theta_v() const;

    //! Degeneracies for each vibrational mode.
    const std::vector<unsigned int>& ndg_v() const;

    //! Characteristic electronic excitation temperatures [K].
    const std::vector<NumericType>& theta_e() const;

    //! Degeneracies for each electronic modes.
    const std::vector<unsigned int>& ndg_e() const;

    //! Formatted print
    /*!
     * Defaults to \p std::cout.
     */
    void print(std::ostream &os = std::cout) const;

    //!Formatted print 
    /*!
     * Allows you to do std::cout << object << std::endl;
     */
    friend std::ostream& operator<<( std::ostream& os,
				     const ChemicalSpecies<NumericType>& species )
    {
      species.print(os);
      return os;
    }

  protected:

    //! Name of chemical species
    const std::string _name;

    //! Molecular weight (or molar mass) in units of [g/mol] or [kg/kmol]
    const NumericType _mol_wght;

    //! Gas constant in units of [J/kg-K]
    const NumericType _R;

    //! Formation enthalpy in units of [J/kg]
    const NumericType _h_form;

    //! Number of translational degrees of freedom
    const NumericType _n_tr_dofs;

    //! Electrical charge number
    const int _charge;

    //! Characteristic vibrational temperature in units of [K]
    std::vector<NumericType> _theta_v;

    //! Degeneracies for each vibrational mode
    std::vector<unsigned int> _ndg_v;

    //! Characteristic electronic temperature in units of [K]
    std::vector<NumericType> _theta_e;

    //! Degeneracies for each electronic mode
    std::vector<unsigned int> _ndg_e;
    
  }; // class ChemicalSpecies


  /* ------------------------- Friend Functions ------------------------- */

  

  /* ------------------------- Inline Functions ------------------------- */

  template<class NumericType>
  inline
  const std::string& ChemicalSpecies<NumericType>::species() const 
  { return _name; }

  template<class NumericType>
  inline
  NumericType ChemicalSpecies<NumericType>::molar_mass() const
  { return _mol_wght; }

  template<class NumericType>
  inline
  NumericType ChemicalSpecies<NumericType>::gas_constant() const 
  { return _R; }

  template<class NumericType>
  inline
  NumericType ChemicalSpecies<NumericType>::formation_enthalpy() const 
  { return _h_form; }

  template<class NumericType>
  inline
  NumericType ChemicalSpecies<NumericType>::n_tr_dofs() const 
  { return _n_tr_dofs; }

  template<class NumericType>
  inline
  int ChemicalSpecies<NumericType>::charge() const 
  { return _charge; }

  template<class NumericType>
  inline
  bool ChemicalSpecies<NumericType>::has_vibrational_modes() const 
  { return !_theta_v.empty(); }

  template<class NumericType>
  inline
  unsigned int ChemicalSpecies<NumericType>::n_vibrational_modes() const
  {
    antioch_assert_equal_to(_theta_v.size(), _ndg_v.size());
    return _theta_v.size();
  }

  template<class NumericType>
  inline
  const std::vector<NumericType>& ChemicalSpecies<NumericType>::theta_v() const
  { return _theta_v; }

  template<class NumericType>
  inline
  const std::vector<unsigned int>& ChemicalSpecies<NumericType>::ndg_v() const 
  { return _ndg_v; }

  template<class NumericType>
  inline
  const std::vector<NumericType>& ChemicalSpecies<NumericType>::theta_e() const 
  { return _theta_e; }

  template<class NumericType>
  inline
  const std::vector<unsigned int>& ChemicalSpecies<NumericType>::ndg_e() const
  { return _ndg_e; }


  template<class NumericType>
  inline
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
  inline
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
  inline
  ChemicalSpecies<NumericType>::~ChemicalSpecies()
  {
    return;
  }


  template<class NumericType>
  inline
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


} // end namespace Antioch

#endif // ANTIOCH_CHEMICAL_SPECIES_H
