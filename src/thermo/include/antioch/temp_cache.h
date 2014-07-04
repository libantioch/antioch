//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014 Paul T. Bauman, Sylvain Plessis, Roy H. Stonger
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

#ifndef ANTIOCH_TEMP_CACHE_H
#define ANTIOCH_TEMP_CACHE_H

namespace Antioch
{
  template<typename StateType=double>
  class TempCache
  {
  public:

    explicit TempCache(const StateType& T_in);

    TempCache(const StateType& T_in, 
              const StateType& T2_in, 
              const StateType& T3_in, 
              const StateType& T4_in, 
              const StateType& lnT_in);

    const StateType& T;
    StateType T2;
    StateType T3;
    StateType T4;
    StateType lnT;

  private:

    TempCache();

  };

  template<typename StateType>
  TempCache<StateType>::TempCache(const StateType& T_in)
    : T(T_in), T2(T*T), T3(T2*T), T4(T2*T2), lnT(T_in)
  {
    using std::log;

    lnT = log(T);
    return;
  }

  template<typename StateType>
  TempCache<StateType>::TempCache(const StateType& T_in, 
                                  const StateType& T2_in, 
                                  const StateType& T3_in, 
                                  const StateType& T4_in, 
                                  const StateType& lnT_in)
    : T(T_in), T2(T2_in), T3(T3_in), T4(T4_in), lnT(lnT_in)
  {
    return;
  }

}

#endif // ANTIOCH_TEMP_CACHE_H
