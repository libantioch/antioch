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
#include "antioch/antioch_asserts.h"

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

        template <typename StateType, typename VIntType>
        StateType custom_bin_value(const StateType & custom_head, const StateType & custom_tail,
                                   const VIntType  & index_heads, unsigned int custom_head_index,
                                   const VectorCoeffType & list_ref_head_tails,
                                   const VectorCoeffType & list_ref_values) const;

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
template <typename VectorStateType>
inline
void SigmaBinConverter<VectorCoeffType>::y_on_custom_grid(const VectorCoeffType &x_old, const VectorCoeffType &y_old,  
                                                          const VectorStateType &x_custom,    VectorStateType &y_custom) const
{
// data consistency
  antioch_assert_not_equal_to(x_custom.size(),0);
  antioch_assert_not_equal_to(x_old.size(),0);
  antioch_assert_equal_to(x_old.size(),y_old.size());

  y_custom.clear();
  y_custom.resize(x_custom.size(),Antioch::zero_clone(x_custom[0]));

// first meta-prog needed stuff

  typedef typename Antioch::value_type<VectorStateType>::type     StateType;
  typedef typename Antioch::rebind<StateType, unsigned int>::type IntType;
  typedef typename Antioch::rebind<VectorStateType,IntType>::type VIntType;

// find all the indexes of old that are just after all
// the custom values
// todo: horrible inefficient way to build the indexes container
  IntType example;
  Antioch::zero_clone(example,x_custom[0]);
  VIntType ihead(x_custom.size(), example);

  for(unsigned int ic = 0; ic < x_custom.size(); ic++)
  {
    IntType ihigh = Antioch::constant_clone(example,x_old.size()-1);
    for (int i = 0; i != x_old.size() - 1; ++i)
    {
      IntType icus  = Antioch::constant_clone(example,ic);
      ihigh = Antioch::if_else (Antioch::eval_index(x_custom,icus) < x_old[i] && ihigh == Antioch::constant_clone(example,x_old.size()-1),
                                Antioch::constant_clone(example,i),
                                ihigh);
    }
    ihead[ic] = ihigh;
  }

  // bin
  for(unsigned int ic = 0; ic < x_custom.size() - 1; ic++) // right stairs, last one = 0
  {
    y_custom[ic] = this->custom_bin_value<StateType,VIntType>(x_custom[ic], x_custom[ic + 1],ihead, ic, x_old, y_old);
  }

  return;

}

   template <typename VectorCoeffType>
   template <typename StateType, typename VIntType>
   inline
   StateType SigmaBinConverter<VectorCoeffType>::custom_bin_value(const StateType & custom_head, const StateType & custom_tail,
                                                                  const VIntType  & index_heads, unsigned int custom_head_index,
                                                                  const VectorCoeffType & list_ref_head_tails,
                                                                  const VectorCoeffType & list_ref_values) const
   {

       using std::min;
       using Antioch::min;
       using std::max;
       using Antioch::max;

       antioch_assert_equal_to(list_ref_head_tails.size(),list_ref_values.size());


       StateType surf = Antioch::zero_clone(custom_head);

        // this is an Antioch::rebind<StateType,unsigned int>::type
       typename Antioch::value_type<VIntType>::type  start_head = index_heads[custom_head_index];
                                                                                // the typename ... ::type> is an unsigned int
       typename Antioch::value_type<VIntType>::type  value_head = Antioch::if_else<typename Antioch::value_type<typename Antioch::value_type<VIntType>::type >::type>
                                                                                (start_head > Antioch::zero_clone(start_head),
                                                                                       start_head - Antioch::constant_clone(start_head,1),
                                                                                       index_heads.back()); // right stairs, value never used

      typename Antioch::value_type<VIntType>::type ref_end_tail = min<typename Antioch::value_type<typename Antioch::value_type<VIntType>::type>::type>
                        (start_head + Antioch::constant_clone(start_head,1),Antioch::constant_clone(start_head,list_ref_head_tails.size() - 1));


        // if StateType is vectorized, it takes care of it here
       StateType ref_head  = Antioch::upgrade_type(custom_head,list_ref_head_tails,start_head);
       StateType ref_tail  = Antioch::upgrade_type(custom_head,list_ref_head_tails,ref_end_tail);
       StateType ref_value = Antioch::upgrade_type(custom_head,list_ref_values,value_head);

       //head from custom head to ref head
       // super not efficient, everything is calculated every time...
       surf += Antioch::if_else(Antioch::constant_clone(custom_head,list_ref_head_tails.front()) > custom_head ||
                                Antioch::constant_clone(custom_head,list_ref_head_tails.back()) < custom_head, // custom is outside ref
                                Antioch::zero_clone(surf),
                                Antioch::if_else<typename Antioch::value_type<StateType>::type>(ref_head < custom_tail,   // custom is within ref bin
                                                 ref_value * (ref_head - custom_head),
                                                 ref_value * (custom_tail - custom_head))
                       );
                                                
       //body from ref head to ref last tail
       while(Antioch::disjunction(ref_tail < custom_tail &&                     // ref is below tail (<=> start_head < index_heads[custom_head_index + 1])
                                  start_head < Antioch::constant_clone(start_head,list_ref_head_tails.size() - 1))) // ref is still defined
       {
           ref_end_tail = min<typename Antioch::value_type<typename Antioch::value_type<VIntType>::type>::type>
                                (start_head + Antioch::constant_clone(start_head,1),Antioch::constant_clone(start_head,list_ref_head_tails.size() - 1));

           ref_head  = Antioch::upgrade_type(custom_head,list_ref_head_tails,start_head);
           ref_tail  = Antioch::upgrade_type(custom_head,list_ref_head_tails,ref_end_tail);
           ref_value = Antioch::upgrade_type(custom_head,list_ref_values,start_head);

           surf += Antioch::if_else<typename Antioch::value_type<StateType>::type>(ref_tail < custom_tail && start_head < Antioch::constant_clone(start_head,list_ref_head_tails.size()),
                                        (ref_tail - ref_head) * ref_value,
                                        Antioch::zero_clone(surf));

           start_head += Antioch::if_else(ref_tail < custom_tail && start_head < Antioch::constant_clone(start_head,list_ref_head_tails.size()),
                                            Antioch::constant_clone(start_head,1),
                                            Antioch::zero_clone(start_head));
       }

      ref_end_tail = min<typename Antioch::value_type<typename Antioch::value_type<VIntType>::type>::type>
                        (start_head + Antioch::constant_clone(start_head,1),Antioch::constant_clone(start_head,list_ref_head_tails.size() - 1));
                                                                                 
       ref_head  = Antioch::upgrade_type(custom_head,list_ref_head_tails,start_head);
       ref_tail  = Antioch::upgrade_type(custom_head,list_ref_head_tails,ref_end_tail);
       ref_value = Antioch::upgrade_type(custom_head,list_ref_values,start_head);

       //tail from ref_head to custom_tail
       surf += Antioch::if_else<typename Antioch::value_type<StateType>::type>(
                        Antioch::constant_clone(custom_tail,list_ref_head_tails.back()) < custom_tail || // custom is outside ref
                        Antioch::constant_clone(custom_tail,list_ref_head_tails.front()) > custom_tail || // custom is outside ref
                                ref_head > custom_tail,   // custom is fully inside ref bin (already taken into account in head)
                                   Antioch::zero_clone(surf),
                                   ref_value * (custom_tail - ref_head));

      return surf / (custom_tail - custom_head);
   }
}

#endif
