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

    /*!
     * @returns species translational specific heat at constant volume, [J/kg-K].
     * Since the translational modes are assumed to be fully polulated
     * this is simply
     * \f[
     *   C^{trans}_{v,s} \equiv \frac{\partial e^{trans}_s}{\partial T} = \frac{3}{2} R_s
     * \f]
     */
    CoeffType cv_trans( const unsigned int species ) const;

    /*!
     * @returns species translational specific over R heat at constant volume.
     * Since the translational modes are assumed to be fully polulated
     * this is simply
     * \f[
     *   \frac{C^{trans}_{v,s}}{\mathrm{R}} = \frac{3}{2}
     * \f]
     */
    CoeffType cv_trans_over_R( const unsigned int species ) const;

    /*!
     * @returns species rotational specific heat at constant volume, [J/kg-K].
     * By convention, we lump the translational/rotational components
     * \f[
     *   C^{tr}_{v,s} \equiv C^{trans}_{v,s} + C^{rot}_{v,s}
     * \f]
     * so then
     * \f[
     *   C^{rot}_{v,s} \equiv C^{tr}_{v,s} - C^{trans}_{v,s}
     * \f]
     */
    CoeffType cv_rot( const unsigned int species ) const;

    /*!
     * @returns species rotational specific heat at constant volume, normalized
     * by the species gas constant.
     * By convention, we lump the translational/rotational components
     * \f[
     *   C^{tr}_{v,s} \equiv C^{trans}_{v,s} + C^{rot}_{v,s}
     * \f]
     * so then
     * \f[
     *   \frac{C^{rot}_{v,s}}{\mathrm{R}} \equiv \frac{C^{tr}_{v,s}}{\mathrm{R}} - \frac{C^{trans}_{v,s}}{\mathrm{R}}
     * \f]
     */
    CoeffType cv_rot_over_R( const unsigned int species ) const;

    /*!
     * @returns species translational+rotational specific heat at
     * constant volume, [J/kg-K].
     */
    CoeffType cv_tr (const unsigned int species) const;

    /*!
     * @returns vector-type of translational+rotational specific heat at
     * constant volume, [J/kg-K], for each corresponding species in the
     * input mass_fractions vector-type.
     */
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type
    >::type
    cv_tr (const VectorStateType& mass_fractions) const;

    /**
     * @returns species translational+rotational specific heat at
     * constant volume, normalized by the species gas constant.
     */
    CoeffType cv_tr_over_R (const unsigned int species) const;

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

  template<typename CoeffType>
  inline
  CoeffType MicroThermoBase<CoeffType>::cv_trans( const unsigned int species ) const
  {
    return CoeffType(1.5)*this->_chem_mixture.R(species);
  }

  template<typename CoeffType>
  inline
  CoeffType MicroThermoBase<CoeffType>::cv_trans_over_R( const unsigned int /*species*/ ) const
  {
    return CoeffType(1.5);
  }

  template<typename CoeffType>
  inline
  CoeffType MicroThermoBase<CoeffType>::cv_rot( const unsigned int species ) const
  {
    using std::max;

    return max(this->cv_tr(species) - this->cv_trans(species), CoeffType(0) );
  }

   template<typename CoeffType>
   inline
  CoeffType MicroThermoBase<CoeffType>::cv_rot_over_R( const unsigned int species ) const
  {
    using std::max;

    return max(this->cv_tr_over_R(species) - this->cv_trans_over_R(species), CoeffType(0) );
  }

  template<typename CoeffType>
  inline
  CoeffType MicroThermoBase<CoeffType>::cv_tr (const unsigned int species) const
  {
    return this->_chem_mixture.R(species)*(this->_chem_mixture.chemical_species()[species])->n_tr_dofs();
  }

  template<typename CoeffType>
  inline
  CoeffType MicroThermoBase<CoeffType>::cv_tr_over_R (const unsigned int species) const
  {
    return (this->_chem_mixture.chemical_species()[species])->n_tr_dofs();
  }

  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type
  >::type
  MicroThermoBase<CoeffType>::cv_tr (const VectorStateType& mass_fractions) const
  {
    typename Antioch::value_type<VectorStateType>::type
      cv_tr = mass_fractions[0]*this->cv_tr(0);

    for( unsigned int s = 1; s < this->_chem_mixture.n_species(); s++ )
      cv_tr += mass_fractions[s]*this->cv_tr(s);

    return cv_tr;
  }

} // end namespace Antioch

#endif // ANTIOCH_MICRO_THERMO_BASE_H
