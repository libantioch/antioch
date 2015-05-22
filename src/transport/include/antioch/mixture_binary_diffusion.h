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

#ifndef ANTIOCH_MIXTURE_BINARY_DIFFUSION_H
#define ANTIOCH_MIXTURE_BINARY_DIFFUSION_H

#include "antioch/diffusion_traits.h"
#include "antioch/mixture_diffusion_base.h"


namespace Antioch
{
  //! Container class for species binary diffusion models
  /*! For the given set of chemical species in the input TransportMixture, this contains
   * the diffusion models for each of those species and provides and interface for
   * computing the species diffusion coefficients. Total diffusion coefficients is computed
   * by a mixing model, e.g. WilkeTransportMixture. This class is templated on the diffusion model,
   * so an inherent assumption is that all species diffusions have the same model.
   */
  template<typename Diffusion, class CoeffType=double>
  class MixtureBinaryDiffusion : public MixtureDiffusionBase<MixtureBinaryDiffusion<Diffusion,CoeffType>,CoeffType>
  {
  public:

    MixtureBinaryDiffusion( const TransportMixture<CoeffType>& transport_mixture );

    virtual ~MixtureBinaryDiffusion();

    //! Add binary diffusion model for species pair i, j
    void add( unsigned int i, unsigned int j,
              const TransportSpecies<CoeffType>& s_i,
              const TransportSpecies<CoeffType>& s_j );

    //! Reset model coefficients for binary diffusion model for species pair i,j
    void reset_coeffs( unsigned int i, unsigned int j,
                       const TransportSpecies<CoeffType>& s_i,
                       const TransportSpecies<CoeffType>& s_j );

    //! Define the diffusion model. Mainly for use in DiffusionTraits
    typedef Diffusion Type;

    //! Friend the base class so we can make the implementation protected
    friend class MixtureDiffusionBase<MixtureBinaryDiffusion<Diffusion,CoeffType>,CoeffType>;

  protected:

    //! Compute the full binary diffusion matrix in D
    /*! \todo We need a matrix type and separate method for exploiting the symmetry*/
    template<typename StateType, typename MatrixStateType>
    void compute_binary_diffusion_matrix_impl( const StateType& T, const StateType& molar_density, MatrixStateType& D ) const;

    //! Should never be called by this class, facilitating interaction in WilkeTransportEvaluator
    /*!
     * There are wildly differing requirements between species and binary diffusion
     * models, but we want to have a uniform interface through WilkeTransportEvaluator.
     * So we need to implement this function, but it will throw an error.
     */
    template<typename StateType>
    void compute_species_diffusivity_impl( unsigned int s, const StateType& rho,
                                           const StateType& cp, const StateType& k,
                                           StateType& D ) const;

    //! Stores binary diffusion model for all species pairs
    /*! \todo This will always be symmetric so we should make a smarter
      container that only stores the upper (lower) triangle of the
      species pairs and use accordingly. */
    std::vector<std::vector<Diffusion*> > _binary_diffusivities;

  private:

    MixtureBinaryDiffusion();
  };

  template<typename Diffusion, class CoeffType>
  MixtureBinaryDiffusion<Diffusion,CoeffType>::MixtureBinaryDiffusion( const TransportMixture<CoeffType>& transport_mixture )
    :  MixtureDiffusionBase<MixtureBinaryDiffusion<Diffusion,CoeffType>,CoeffType>(transport_mixture),
    _binary_diffusivities( transport_mixture.n_species() )
  {
#ifdef ANTIOCH_HAVE_CXX_STATIC_ASSERT
    static_assert( DiffusionTraits<Diffusion,CoeffType>::is_binary_diffusion,
                   "Can only instantiate MixtureBinaryDiffusion with a binary diffusion model!" );
#endif

    // Finish allocating space for binary diffusion coeffient objects
    for( unsigned int s = 0; s < transport_mixture.n_species(); s++ )
      _binary_diffusivities[s].resize( transport_mixture.n_species(), NULL );
  }

  template<typename Diffusion, class CoeffType>
  MixtureBinaryDiffusion<Diffusion,CoeffType>::~MixtureBinaryDiffusion()
  {
    for( typename std::vector<std::vector<Diffusion*> >::iterator it_outer = _binary_diffusivities.begin();
         it_outer != _binary_diffusivities.end(); ++it_outer )
      for( typename std::vector<Diffusion*>::iterator it_inner = it_outer->begin();
           it_inner != it_outer->end(); ++it_inner )
        delete *it_inner;
  }

  template<typename Diffusion, class CoeffType>
  void MixtureBinaryDiffusion<Diffusion,CoeffType>::add( unsigned int i, unsigned int j,
                                                         const TransportSpecies<CoeffType>& s_i,
                                                         const TransportSpecies<CoeffType>& s_j )
  {
    antioch_assert_less( i, _binary_diffusivities.size() );
    antioch_assert_less( j, _binary_diffusivities[i].size() );

    _binary_diffusivities[i][j] = new Diffusion( s_i, s_j );
  }

  template<typename Diffusion, class CoeffType>
  void MixtureBinaryDiffusion<Diffusion,CoeffType>::reset_coeffs( unsigned int i, unsigned int j,
                                                                  const TransportSpecies<CoeffType>& s_i,
                                                                  const TransportSpecies<CoeffType>& s_j )
  {
    antioch_assert_less( i, _binary_diffusivities.size() );
    antioch_assert_less( j, _binary_diffusivities[i].size() );
    antioch_assert(_binary_diffusivities[i][j]);

    _binary_diffusivities[i][j]->reset_coeffs( s_i, s_j );
  }

  template<typename Diffusion, class CoeffType>
  template<typename StateType, typename MatrixStateType>
  inline
  void MixtureBinaryDiffusion<Diffusion,CoeffType>::compute_binary_diffusion_matrix_impl( const StateType& T,
                                                                                          const StateType& molar_density,
                                                                                          MatrixStateType& D ) const
  {
    const unsigned int n_cols = D.size();
    antioch_assert_greater(n_cols,0);

    // Make sure it's a square matrix
#ifndef NDEBUG
    for( unsigned int r = 0; r < n_cols; r++ )
      antioch_assert_equal_to(D[r].size(),n_cols);
#endif

    for(unsigned int i = 0; i < D.size(); i++)
      {
        for(unsigned int j = 0; j < D.size(); j++)
          {
            antioch_assert(_binary_diffusivities[i][j]);

            D[i][j] = (*_binary_diffusivities[i][j])(T,molar_density);
          }
      }
  }

  template<typename Diffusion, class CoeffType>
  template<typename StateType>
  void MixtureBinaryDiffusion<Diffusion,CoeffType>::compute_species_diffusivity_impl( unsigned int s,
                                                                                      const StateType& rho,
                                                                                      const StateType& cp,
                                                                                      const StateType& k,
                                                                                      StateType& D ) const
  {
    std::string error = "ERROR: You're trying to use a species diffusion implementation\n";
    error += "       with a binary diffusion model!\n";
    antioch_msg_error(error);
  }

} // end namespace Antioch

#endif // ANTIOCH_MIXTURE_BINARY_DIFFUSION_H
