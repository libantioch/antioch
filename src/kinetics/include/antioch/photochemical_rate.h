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

#ifndef _PHOTOCHEMICAL_RATE_
#define _PHOTOCHEMICAL_RATE_

//C++
#include <vector>

//Antioch
#include "antioch/cmath_shims.h"
#include "antioch/metaprogramming.h"

namespace Antioch{

template<typename StateType, typename VectorStateType>
class PhotoRate{
     public:
       PhotoRate(const VectorStateType &cs, const VectorStateType &lambda):
                _cross_section(cs),_lambda_grid(lambda){}
       ~PhotoRate(){}
       StateType forward_rate_constant(const VectorStateType &hv, const VectorStateType &lambda);

       void set_cross_section(const VectorStateType &cs) {_cross_section = cs;}
       void set_lambda_grid(const VectorStateType &l)    {_lambda_grid = l;}

     private:
       VectorStateType _cross_section;
       VectorStateType _lambda_grid;

       VectorStateType make_grid_cool(const VectorStateType &hv, const VectorStateType &lambda);

};

template<typename StateType, typename VectorStateType>
StateType PhotoRate<StateType,VectorStateType>::forward_rate_constant(const VectorStateType &hv, const VectorStateType &lambda)
{
// lambda grid
  VectorStateType hvgrid = this->make_grid_cool(hv,lambda);
// integration, those are bins => just multiply
  StateType rfwd;
  Antioch::set_zero(rfwd);
  for(unsigned int i = 0; i < hvgrid.size(); i++)
  {
    rfwd += _cross_section[i] * hvgrid[i];
  }

  return rfwd;
}

template <typename StateType, typename VectorStateType>
VectorStateType PhotoRate<StateType,VectorStateType>::make_grid_cool(const VectorStateType &hv, const VectorStateType &lambda)
{
  VectorStateType out_hv_on_grid;
  out_hv_on_grid.resize(_lambda_grid.size(),0.L);
  unsigned int j(1);

  for(unsigned int i = 0; i < _lambda_grid.size()-1; i++)//bin per bin, bin hv[j] between lambda[j] and lambda[j+1]
  {
     while(lambda[j] < _lambda_grid[i])
     {
       j++;
       if(j >= lambda.size())return out_hv_on_grid;
     }
     StateType bin;
     Antioch::set_zero(bin);
     StateType diff_min = lambda[j] - _lambda_grid[i];
     StateType diff_max = _lambda_grid[i+1] - lambda[j];
     bin += hv[j]   * (diff_min)/(lambda[j] - lambda[j-1]);
     if(j < lambda.size() - 1)bin += hv[j+1] * (diff_max)/(lambda[j+1] - lambda[j]);
     out_hv_on_grid[i] = bin;
  }

  return out_hv_on_grid;
}

}

#endif
