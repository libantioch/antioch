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

#ifndef _ANTIOCH_LINDEMANN_FALLOFF_H
#define _ANTIOCH_LINDEMANN_FALLOFF_H

//Antioch

//C++

namespace Antioch{
/*!\class LindemannFalloff
 *
 * The Lindemann model is the simplest falloff
 * model:
 * \f[
 *     F = 1
 * \f]
 */
template <typename CoeffType = double>
class LindemannFalloff{
     public:
       LindemannFalloff();
       ~LindemannFalloff();

     CoeffType compute_F(const CoeffType& T, const CoeffType &Pr) const;

};
  template<typename CoeffType>
  inline
  CoeffType LindemannFallOff<CoeffType>::compute_F(const CoeffType &T, const CoeffType &Pr) const
  {
    return 1.;
  }

  template<typename CoeffType>
  inline
  LindemannFallOff<CoeffType>::LindemannFalloff()
  {
    return;
  }

  template<typename CoeffType>
  inline
  LindemannFallOff<CoeffType>::~LindemannFalloff()
  {
    return;
  }
}
