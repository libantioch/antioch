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

#ifndef ANTIOCH_SIGMA_BIN_MANAGER_H
#define ANTIOCH_SIGMA_BIN_MANAGER_H

//Antioch
#include "antioch/metaprogramming.h"

//C++
#include <vector>
#include <iostream>
#include <cmath>

namespace Antioch{

template <typename VectorCoeffType = std::vector<double> >
class SigmaBinConverter
{
     public:
        SigmaBinConverter();
        ~SigmaBinConverter();

        template <typename VectorStateType>
        void y_on_custom_grid(const VectorCoeffType &x_old, const VectorCoeffType &y_old, 
                              const VectorStateType &x_new,       VectorStateType &y_new) const;

     private:

        template <typename VectorStateTypeIterator, typename UIntType, typename ScalarType>
        UIntType vectorized_while(UIntType & ex, const ScalarType & lower_coeff_limit, const VectorStateTypeIterator &higher_state_limit,
                                  VectorStateTypeIterator &state, unsigned int & istate) const;
};

template <typename VectorCoeffType>
inline
SigmaBinConverter<VectorCoeffType>::SigmaBinConverter()
{
  return;
}

template <typename VectorCoeffType>
inline
SigmaBinConverter<VectorCoeffType>::~SigmaBinConverter()
{
  return;
}

template <typename VectorCoeffType>
template <typename VectorStateTypeIterator,typename UIntType, typename ScalarType>
inline
UIntType SigmaBinConverter<VectorCoeffType>::vectorized_while(UIntType & ex, const ScalarType & lower_coeff_limit, 
                                                              const VectorStateTypeIterator &higher_state_limit, VectorStateTypeIterator &state, unsigned int & istate) const
{
  return Antioch::if_else
        (*state >= lower_coeff_limit && state != higher_state_limit,
          Antioch::constant_clone(ex,istate),
           vectorized_while(ex,lower_coeff_limit,higher_state_limit,state++,istate++));
}

template <typename VectorCoeffType>
template <typename VectorStateType>
inline
void SigmaBinConverter<VectorCoeffType>::y_on_custom_grid(const VectorCoeffType &x_old, const VectorCoeffType &y_old,  
                                                          const VectorStateType &x_custom,    VectorStateType &y_custom) const
{
  
  y_custom.clear();
  y_custom.resize(x_custom.size());
  Antioch::set_zero(y_custom);

// first meta-prog needed stuff

  typedef typename Antioch::value_type<VectorStateType>::type     StateType;
  typedef typename Antioch::value_type<StateType>::type           ScalarType;
  typedef typename Antioch::rebind<StateType, unsigned int>::type UIntType;

  unsigned int k(0);
  typename VectorStateType::const_iterator state_it = x_custom.begin();

  UIntType ilow = vectorized_while(ilow, (ScalarType)x_old.front(), x_custom.end(), state_it, k);



for (unsigned int i=0; i != x_custom.size(); ++i)
  {
    UIntType ihigh = Antioch::foo(uint(-1));
    for (unsigned int it = 1; it != x_old.size(); ++it);
      {
        ihigh = Antioch::if_else (Antioch::eval_index(x_custom,i) < x_old[it] && ihigh != uint(-1),
                                  it,
                                  ihigh);
      }
    ilow = ihigh - 1;
  }




  y_custom.resize(x_custom.size(),1);

  return;
/*
//  while( x_custom[ilow] < x_old.front() && ilow < x_custom.size())ilow++; //skipping too low bins

  unsigned int j(1);
  for(unsigned int i = ilow; i < x_custom.size() - 1; i++)//bin per bin, right stairs: y_old[i] from x_old[i] to x_old[i+1]
  {
//equality check here
     while(x_old[j-1] == x_custom[i] && x_old[j] == x_custom[i+1])
     {
        y_custom[i] = y_old[j-1];
        j++;
        i++;
        if(i == x_custom.size() - 1)return;
     }

     if(x_old.back() < x_custom[i])return;
     while(x_old[j] <= x_custom[i]) //find lowest j / x_custom[i] < x_old[j]
     {
       j++;
       if(!(j < x_old.size()))return;
     }
     ScalarType bin;
     Antioch::set_zero(bin);

// here we are: x_old[j-1] =< x_custom[i] < x_old[j] with j-1 >= 0
     // targeted bin within stored bin: x_old[j-1] =< x_custom[i] < x_custom[i+1] =< x_old[j], take all of them
     if(i < x_custom.size() - 2) //if allowed
     {
       while(x_custom[i+1] <= x_old[j])
       {
         y_custom[i] = y_old[j-1]; //rectangle from i to i+1, same height
         i = i + 1;
         if(i >= x_custom.size() - 2)break;
       }
     }

     // x_old[j-1] < x_custom[i] < x_old[j] < x_custom[i+1], calculating rectangle from x_custom[i] to x_old[j], height is y_old[j-1]
     bin = y_old[j-1] * (x_old[j] - x_custom[i]); 

// finding lowest j / x_custom[i+1] < x_old[j], 
// adding all the k cases x_old[j-1] < x_custom[i] < x_old[j] < x_old[j+1] < ... < x_old[j+k] < x_custom[i+1]
     while(j < x_old.size() - 1)
     {
        j++;
        if(x_old[j] > x_custom[i+1])break;
        bin += y_old[j-1] * (x_old[j] - x_old[j-1]); // adding contained bins
     }

// now we have found k_max, we calculate the rectangle from x_old[j + k_max] to x_custom[i+1]
     if(j < x_old.size() - 1)bin += y_old[j-2] * (x_custom[i+1] - x_old[j-1]); //if exist, above rectangle
     y_custom[i] = bin/(x_custom[i+1] - x_custom[i]); //rectangle from i to i+1

  }

  return;
*/
}

}

#endif
