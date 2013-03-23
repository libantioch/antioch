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

#ifndef ANTIOCH_CONSTANT_LEWIS_DIFFUSIVITY_H
#define ANTIOCH_CONSTANT_LEWIS_DIFFUSIVITY_H

namespace Antioch
{

  template<typename CoeffType=double>
  class ConstantLewisDiffusivity
  {
  public:
    
    ConstantLewisDiffusivity( const CoeffType Le );

    ~ConstantLewisDiffusivity();
    
    template<typename StateType>
    StateType D( const StateType rho, const StateType cp, const StateType k ) const;

  protected:

    CoeffType _Le;
    
  };

  template<typename CoeffType>
  ConstantLewisDiffusivity<CoeffType>::ConstantLewisDiffusivity( const CoeffType Le )
    : _Le(Le)
  {
    return;
  }

  template<typename CoeffType>
  ConstantLewisDiffusivity<CoeffType>::~ConstantLewisDiffusivity()
  {
    return;
  }

  template<typename CoeffType>
  template<typename StateType>
  StateType ConstantLewisDiffusivity<CoeffType>::D( const StateType rho,
                                                    const StateType cp,
                                                    const StateType k ) const
  {
    return _Le*k/(rho*cp);
  }

} // end namespace Antioch

#endif // ANTIOCH_CONSTANT_LEWIS_DIFFUSIVITY_H
