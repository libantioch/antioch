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

template<typename CoeffType, typename VectorCoeffType>
class PhotoRate{

     private:
       VectorCoeffType _cross_section;
       VectorCoeffType _lambda_grid;

       VectorCoeffType harmonize_grid(const VectorCoeffType &hv, const VectorCoeffType &lambda);

     public:
       PhotoRate(const VectorCoeffType &cs, const VectorCoeffType &lambda):
                _cross_section(cs),_lambda_grid(lambda){}
       ~PhotoRate(){}


       CoeffType forward_rate_constant(const VectorCoeffType &hv, const VectorCoeffType &lambda);

       void set_cross_section(const VectorCoeffType &cs) {_cross_section = cs;}
       void set_lambda_grid(const VectorCoeffType &l)    {_lambda_grid = l;}

};

template<typename CoeffType, typename VectorCoeffType>
CoeffType PhotoRate<CoeffType,VectorCoeffType>::forward_rate_constant(const VectorCoeffType &hv, const VectorCoeffType &lambda)
{
// lambda grid
  VectorCoeffType hvgrid = this->harmonize_grid(hv,lambda);
// integration, those are bins => just multiply
  CoeffType rfwd;
  Antioch::set_zero(rfwd);
  for(unsigned int i = 0; i < hvgrid.size(); i++)
  {
    rfwd += _cross_section[i] * hvgrid[i];
  }

  return rfwd;
}

template <typename CoeffType, typename VectorCoeffType>
VectorCoeffType PhotoRate<CoeffType,VectorCoeffType>::harmonize_grid(const VectorCoeffType &hv, const VectorCoeffType &lambda)
{
  VectorCoeffType out_hv_on_grid;
  out_hv_on_grid.resize(_lambda_grid.size(),0.L);
  unsigned int j(1);

  for(unsigned int i = 0; i < _lambda_grid.size()-1; i++)//bin per bin, bin hv[j] between lambda[j] and lambda[j+1]
  {
     while(lambda[j] < _lambda_grid[i])
     {
       j++;
       if(j >= lambda.size())return out_hv_on_grid;
     }
     CoeffType bin;
     Antioch::set_zero(bin);
     CoeffType diff_min = lambda[j] - _lambda_grid[i];
     CoeffType diff_max = _lambda_grid[i+1] - lambda[j];
     bin += hv[j]   * (diff_min)/(lambda[j] - lambda[j-1]);
     if(j < lambda.size() - 1)bin += hv[j+1] * (diff_max)/(lambda[j+1] - lambda[j]);
     out_hv_on_grid[i] = bin;
  }

  return out_hv_on_grid;
}

}

#endif
