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

#ifndef ANTIOCH_CHEMICAL_MIXTURE_H
#define ANTIOCH_CHEMICAL_MIXTURE_H

// Antioch
#include "antioch_config.h" // default files
#include "antioch/antioch_asserts.h"
#include "antioch/chemical_species.h"
#include "antioch/default_filename.h"
#include "antioch/metaprogramming.h"
#include "antioch/species_parsing.h"
#include "antioch/ascii_parser.h"

// C++
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <iomanip>

namespace Antioch
{
  typedef unsigned int Species;

  //! Class storing chemical mixture properties
  /*!
    This class manages the list of ChemicalSpecies for a requested set
    of species from input.

    We require four types of data:
       - a list of species, in a form of std::vector<std::string>
       - a file for mandatory species information:
                        * molar mass
                        * enthalpie of formation @ 0 K
                        * number of translational-rotational degrees of freedom
                        * charge
       - a file for optional vibrational data
                        * vibrational levels and associated degenerescence
       - a file for optional electronic data
                        * electronic levels and associated degenerescence
    Constructor is overloaded so all or any data can be given, in the
    order above.
  */
  template<typename CoeffType=double>
  class ChemicalMixture
  {
  public:
    //! ascii parser by default
    ChemicalMixture( const std::string & filename = DefaultFilename::species_list(),
                     const bool verbose = true,
                     const std::string & species_data = DefaultFilename::chemical_mixture(),
                     const std::string & vibration_data = DefaultFilename::vibrational_data(),
                     const std::string & electronic_data = DefaultFilename::electronic_data());

    //! ascii parser by default
    ChemicalMixture( const std::vector<std::string>& species_list,
                     const bool verbose = true,
                     const std::string & species_data = DefaultFilename::chemical_mixture(),
                     const std::string & vibration_data = DefaultFilename::vibrational_data(),
                     const std::string & electronic_data = DefaultFilename::electronic_data());

    //! explicit parser, file in parser must contains species list
    ChemicalMixture( ParserBase<CoeffType> * parser,
                     const std::string & species_data = DefaultFilename::chemical_mixture(),
                     const std::string & vibration_data = DefaultFilename::vibrational_data(),
                     const std::string & electronic_data = DefaultFilename::electronic_data());

    ~ChemicalMixture();

    /*! method to send this back

        added for backward compatibility for WilkeMixture v/s
        MixtureAveragedTransportMixture, the idea is to send a chemical mixture
        whatever mixture is provided
     */
    const ChemicalMixture<CoeffType> & chemical_mixture() const;

    //! method to initialize, backward compatibility
    void initialize_species(const std::vector<std::string> & species_list);

    //! method to read characteristics, using one parser
    void read_species_characteristics(ParserBase<CoeffType> * parser,
                                      const std::string & species_data,
                                      const std::string & vibration_data,
                                      const std::string & electronic_data);

    //! method to read mandatory characteristics
    void read_species_mandatory_characteristics(ParserBase<CoeffType> * parser);

    //! method to read vibrational characteristics
    void read_species_vibrational_characteristics(ParserBase<CoeffType> * parser);

    //! method to read electronic characteristics
    void read_species_electronic_characteristics(ParserBase<CoeffType> * parser);

    //! Returns the number of species in this mixture.
    unsigned int n_species() const;

    void add_species( const unsigned int index, const std::string& name,
		      CoeffType mol_wght, CoeffType h_form,
		      CoeffType n_tr_dofs, int charge );

    void add_species_vibrational_data( const unsigned int index,
                                       const CoeffType theta_v,
                                       const unsigned int ndg_v );

    void add_species_electronic_data( const unsigned int index,
                                      const CoeffType theta_e,
                                      const unsigned int ndg_e );

    const std::vector<ChemicalSpecies<CoeffType>*>& chemical_species() const;

    const std::vector<Species>& species_list() const;

    const std::map<std::string,Species>& species_name_map() const;

    const std::map<Species,std::string>& species_inverse_name_map() const;

    const std::map<std::string,Species>& active_species_name_map() const;

    //! Gas constant for species s in [J/kg-K]
    CoeffType R( const unsigned int s ) const;

    //! Gas constant for mixture in [J/kg-K]
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type
    >::type
    R( const VectorStateType& mass_fractions ) const;

    //! Molecular weight (molar mass) for species s in kg/mol
    CoeffType M( const unsigned int s ) const;

    //! Molecular weight (molar mass) for mixture in kg/mol
    /*!
      \f$ \frac{1}{M} = \sum_s \frac{w_s}{M_s}\f$ where
      \f$ w_s \f$ is the mass fraction of species \f$ s \f$ and
      \f$ M_s \f$ is the molecular weight (molar mass) of species \f$ s \f$
    */
    template<typename VectorStateType>
    typename enable_if_c<
      has_size<VectorStateType>::value,
      typename Antioch::value_type<VectorStateType>::type
    >::type
    M( const VectorStateType& mass_fractions ) const;

    //! Species mole fraction
    /*!
      Given mixture molar mass M and mass fraction for species,
      compute species mole fraction using the relationship
      \f$ w_i = x_i \frac{M_i}{M} \f$
    */
    template<typename StateType>
    ANTIOCH_AUTO(StateType)
    X( const unsigned int species, const StateType& M,
       const StateType& mass_fraction ) const
    ANTIOCH_AUTOFUNC(StateType, mass_fraction*M/this->M(species))

    //! All species mole fractions
    /*!
      The output argument mole_fractions should already be properly
      sized to hold the output.
    */
    template<typename VectorStateType>
    void X( typename Antioch::value_type<VectorStateType>::type M,
	    const VectorStateType& mass_fractions,
	    VectorStateType& mole_fractions ) const;

    //! Species molar density
    /*!
      Given total density rho and mass fraction for species,
      compute species moles per unit volume
    */
    template<typename StateType>
    ANTIOCH_AUTO(StateType)
    molar_density( const unsigned int species,
		   const StateType& rho,
		   const StateType& mass_fraction ) const
    ANTIOCH_AUTOFUNC(StateType, rho*mass_fraction/this->M(species))

    //! Species molar densities
    /*!
      Given total density rho and mass fractions for all species,
      compute moles per unit volume for all species
    */
    template<typename StateType, typename VectorStateType>
    void molar_densities( const StateType& rho,
			  const VectorStateType& mass_fractions,
			  VectorStateType& molar_densities ) const;

  protected:

    void init_species_name_map(const std::vector<std::string> & species_list);
    void build_inverse_name_map();

    std::vector<Species> _species_list;
    std::vector<ChemicalSpecies<CoeffType>*> _chemical_species;
    std::map<std::string,Species> _species_name_map;
    std::map<Species,std::string> _species_inv_name_map;

  };


  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType>
  inline
  unsigned int ChemicalMixture<CoeffType>::n_species() const
  {
    return _species_list.size();
  }

  template<typename CoeffType>
  inline
  const std::vector<Species>& ChemicalMixture<CoeffType>::species_list() const
  {
    return _species_list;
  }

  template<typename CoeffType>
  inline
  const std::vector<ChemicalSpecies<CoeffType>*>& ChemicalMixture<CoeffType>::chemical_species() const
  {
    return _chemical_species;
  }

  template<typename CoeffType>
  inline
  const std::map<std::string,Species>& ChemicalMixture<CoeffType>::species_name_map() const
  {
    return _species_name_map;
  }

  template<typename CoeffType>
  inline
  const std::map<Species,std::string>& ChemicalMixture<CoeffType>::species_inverse_name_map() const
  {
    return _species_inv_name_map;
  }

  template<typename CoeffType>
  inline
  const std::map<std::string,Species>&
  ChemicalMixture<CoeffType>::active_species_name_map() const
  {
    antioch_deprecated();
    return _species_name_map;
  }

  template<typename CoeffType>
  inline
  CoeffType ChemicalMixture<CoeffType>::R( const unsigned int s ) const
  {
    return (_chemical_species[s])->gas_constant();
  }

  template<typename CoeffType>
  inline
  CoeffType ChemicalMixture<CoeffType>::M( const unsigned int s ) const
  {
    return (_chemical_species[s])->molar_mass();
  }

  template<typename CoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void ChemicalMixture<CoeffType>::molar_densities( const StateType& rho,
						    const VectorStateType& mass_fractions,
						    VectorStateType& molar_densities ) const
  {
    antioch_assert_equal_to( mass_fractions.size(), this->n_species() );
    antioch_assert_equal_to( molar_densities.size(), this->n_species() );

    // !\todo Figure out how to make this assert vector-compatible
    // antioch_assert_greater( rho, 0.0 );

    for( unsigned int s = 0; s < this->n_species(); s++ )
      {
	molar_densities[s] = rho*mass_fractions[s]/this->M(s);
      }
    return;
  }

  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type
  >::type
  ChemicalMixture<CoeffType>::R( const VectorStateType& mass_fractions ) const
  {
    antioch_assert_equal_to( mass_fractions.size(), _chemical_species.size() );
    antioch_assert_greater( mass_fractions.size(), 0);

    typename Antioch::value_type<VectorStateType>::type
      R = mass_fractions[0]*this->R(0);
    for( unsigned int s = 1; s < mass_fractions.size(); s++ )
      {
	R += mass_fractions[s]*this->R(s);
      }

    return R;
  }


  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  typename enable_if_c<
    has_size<VectorStateType>::value,
    typename Antioch::value_type<VectorStateType>::type
  >::type
  ChemicalMixture<CoeffType>::M( const VectorStateType& mass_fractions ) const
  {
    antioch_assert_equal_to( mass_fractions.size(), _chemical_species.size() );

    typename Antioch::value_type<VectorStateType>::type
      M = mass_fractions[0]/this->M(0);
    for( unsigned int s = 1; s < mass_fractions.size(); s++ )
      {
	M += mass_fractions[s]/(this->M(s));
      }

    return CoeffType(1)/M;
  }


  template<typename CoeffType>
  template<typename VectorStateType>
  inline
  void ChemicalMixture<CoeffType>::X( typename Antioch::value_type<VectorStateType>::type M,
				      const VectorStateType& mass_fractions,
				      VectorStateType& mole_fractions ) const
  {
    antioch_assert_equal_to( mass_fractions.size(), _chemical_species.size() );

    // Require output size to be correct.  This is more efficient and
    // works around the difficulty of resizing Eigen containers of
    // types which may lack default constructors.
    antioch_assert_equal_to( mass_fractions.size(), mole_fractions.size() );

    for( unsigned int s = 0; s < mass_fractions.size(); s++ )
      {
	mole_fractions[s] = this->X(s, M, mass_fractions[s]);
      }

    return;
  }

  template<typename CoeffType>
  inline
  const ChemicalMixture<CoeffType> & ChemicalMixture<CoeffType>::chemical_mixture() const
  {
     antioch_deprecated();
     return *this;
  }

} // end namespace Antioch

#endif //ANTIOCH_CHEMICAL_MIXTURE_H
