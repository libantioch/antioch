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

#ifndef ANTIOCH_LENNARD_JONES_POTENTIAL_H
#define ANTIOCH_LENNARD_JONES_POTENTIAL_H

//Antioch

//C++

namespace Antioch{

  /**
   * The Lennard-Jones potential, for the moment, simply
   * a storage faiclity, easily expandable to a full
   * analysis object
   * \f[
   *    LJ(r_{ij}) = 4\epsilon_{ij} \left[\left(\frac{\sigma_{ij}}{r}\right)^{12} - \left(\frac{\sigma_{ij}}{r}\right)^6\right]
   * \f]
   * with \f$\epsilon\f$ the well depth and \f$\sigma\f$ the collision diameter.
   * The pair values are related to species values by the relations
   * \f[
   *     \epsilon_{ij} = \sqrt{\epsilon_i \epsilon_j}
   * \f]
   * and
   * \f[
   *     \sigma_{ij} = \frac{\sigma_i + \sigma_j}{2}
   * \f]
   *
   */
  template <typename CoeffType>
  class LennardJonesPotential
  {
      public:

        LennardJonesPotential(const CoeffType & depth = 0., const CoeffType & diameter = 0.);
        ~LennardJonesPotential();

        template <typename StateType>
        void set_diameter(const StateType & diameter);

        template <typename StateType>
        void set_depth(const StateType & depth);

        template <typename StateType>
        void reset_coeffs(const StateType & depth, const StateType & diameter);

        CoeffType diameter() const;

        CoeffType depth() const;

      private:


        CoeffType _depth;
        CoeffType _diameter;

  };

  template <typename CoeffType>
  inline
  LennardJonesPotential<CoeffType>::LennardJonesPotential(const CoeffType & depth, const CoeffType & diameter):
        _depth(depth),
        _diameter(diameter)
  {
      return;
  }

  template <typename CoeffType>
  inline
  LennardJonesPotential<CoeffType>::~LennardJonesPotential()
  {
     return;
  }

  template <typename CoeffType>
  inline
  CoeffType LennardJonesPotential<CoeffType>::depth() const
  {
    return _depth;
  }

  template <typename CoeffType>
  inline
  CoeffType LennardJonesPotential<CoeffType>::diameter() const
  {
    return _diameter;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void LennardJonesPotential<CoeffType>::set_diameter(const StateType & diameter)
  {
     _diameter = diameter;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void LennardJonesPotential<CoeffType>::set_depth(const StateType & depth)
  {
     _depth = depth;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void LennardJonesPotential<CoeffType>::reset_coeffs(const StateType & depth, const StateType & diameter)
  {
      this->set_depth(depth);
      this->set_diameter(diameter);
  }
  

} //end namespace Antioch

#endif
