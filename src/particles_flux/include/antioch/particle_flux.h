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

#ifndef ANTIOCH_PARTICLE_FLUX_H
#define ANTIOCH_PARTICLE_FLUX_H

namespace Antioch
{

  /*!\class ParticleFlux
   * Stores the incoming flux of particles
   *
   */
  template<typename VectorCoeffType>
  class ParticleFlux
  {
     private:
        VectorCoeffType _abscissa;
        VectorCoeffType _flux;
        bool _updated;
        bool _x_updated;
        unsigned int _n_coupled;
        unsigned int _n_updated;

     public:
        ParticleFlux();
        ParticleFlux(const VectorCoeffType &x, const VectorCoeffType &flux);
        ~ParticleFlux();

        //!
        bool updated() const;

        //!
        bool x_updated() const;

        //!
        const VectorCoeffType &abscissa() const;

        //!
        const VectorCoeffType &flux() const;

        //!
        template<typename VectorStateType>
        void set_abscissa(const VectorStateType &x);

        //!
        template<typename VectorStateType>
        void set_flux(const VectorStateType &flux);

        //!
        void update_done();

        //!
        void add_a_reaction();

        //!
        void remove_a_reaction();

  };

  template<typename VectorCoeffType>
  inline
  const VectorCoeffType &ParticleFlux<VectorCoeffType>::abscissa() const
  {
     return _abscissa;
  }

  template<typename VectorCoeffType>
  inline
  const VectorCoeffType &ParticleFlux<VectorCoeffType>::flux() const
  {
     return _flux;
  }

  template<typename VectorCoeffType>
  inline
  void ParticleFlux<VectorCoeffType>::remove_a_reaction()
  {
     _n_coupled--;
     return;
  }

  template<typename VectorCoeffType>
  inline
  void ParticleFlux<VectorCoeffType>::add_a_reaction()
  {
     _n_coupled++;
     _updated = true;
     return;
  }

  template<typename VectorCoeffType>
  inline
  bool ParticleFlux<VectorCoeffType>::updated() const
  {
     return _updated;
  }

  template<typename VectorCoeffType>
  inline
  bool ParticleFlux<VectorCoeffType>::x_updated() const
  {
     return _x_updated;
  }

  template<typename VectorCoeffType>
  inline
  void ParticleFlux<VectorCoeffType>::update_done()
  {
     _n_updated++;
     if(_n_updated == _n_coupled)
     {
        _updated   = false;
        _x_updated = false;
        _n_updated = 0;
     }
     return;
  }

  template<typename VectorCoeffType>
  template<typename VectorStateType>
  inline
  void ParticleFlux<VectorCoeffType>::set_abscissa(const VectorStateType &x)
  {
     _abscissa = x;
     _updated = true;
     _x_updated = true;
  }

  template<typename VectorCoeffType>
  template<typename VectorStateType>
  inline
  void ParticleFlux<VectorCoeffType>::set_flux(const VectorStateType &flux)
  {
     _flux = flux;
     _updated = true;
  }

  template<typename VectorCoeffType>
  inline
  ParticleFlux<VectorCoeffType>::ParticleFlux():
  _updated(false),
  _x_updated(false),
  _n_coupled(0),
  _n_updated(0)
  {
    return;
  }

  template<typename VectorCoeffType>
  inline
  ParticleFlux<VectorCoeffType>::~ParticleFlux()
  {
    return;
  }

  template<typename VectorCoeffType>
  inline
  ParticleFlux<VectorCoeffType>::ParticleFlux(const VectorCoeffType &x, const VectorCoeffType &flux):
  _abscissa(x),
  _flux(flux),
  _updated(true),
  _x_updated(true),
  _n_coupled(0),
  _n_updated(0)
  {
    return;
  }

}

#endif
