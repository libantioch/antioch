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

#ifndef ANTIOCH_DATA_EQUILIBRIUM_H
#define ANTIOCH_DATA_EQUILIBRIUM_H

// antioch
#include "antioch/chemical_mixture.h"

//C++
#include <vector>

namespace Antioch
{

//forward declaration
  template<typename CoeffType>
  class ReactionSet;

  template<typename CoeffType = double>
  class DataEquilibrium
  {
     public:
          DataEquilibrium(const CoeffType &T_mix, const CoeffType &P_mix, 
                    const ReactionSet<CoeffType> &reac_set);
          virtual ~DataEquilibrium();

          virtual void fill_constrain(const std::vector<CoeffType> &molar_densities,
                                        std::vector<CoeffType> &F,std::vector<std::vector<CoeffType> > &jacob){return;}// no constrain here

          const CoeffType T() const {return _T;}
          const CoeffType P() const {return _P;}
          virtual const CoeffType local_pressure(const std::vector<CoeffType> &molar_densities) const; //no contrain on P, ideal gas
          const ReactionSet<CoeffType> reaction_set() const {return _reac_set;}

          virtual const unsigned int n_constrain() const {return 0;}
          virtual void update_constrain(unsigned int icstr, const CoeffType &delta_cstr) {return;} //no constrain


     protected:
     const CoeffType &_T;
     const CoeffType &_P;
     const ReactionSet<CoeffType> &_reac_set;

     private:
//I don't want it to be used, change it if you want
     DataEquilibrium(){return;}
  };


  template<typename CoeffType>
  inline
  DataEquilibrium<CoeffType>::DataEquilibrium(const CoeffType &T_mix, const CoeffType &P_mix, 
                    const ReactionSet<CoeffType> &reac_set):
     _T(T_mix),_P(P_mix),_reac_set(reac_set)
  {
     return;
  }

  template<typename CoeffType>
  inline
  DataEquilibrium<CoeffType>::~DataEquilibrium()
  {
    return;
  }

  template<typename CoeffType>
  inline
  const CoeffType DataEquilibrium<CoeffType>::local_pressure(const std::vector<CoeffType> &molar_densities) const
  {
     antioch_assert_equal_to(molar_densities.size(),this->_reac_set.n_species());
     CoeffType sum_mol;
     Antioch::set_zero(sum_mol);
     for(unsigned int i = 0; i < this->_reac_set.n_species(); i++)
     {
        sum_mol += molar_densities[i];
     }
     return (sum_mol * this->_T * Constants::R_universal<CoeffType>());
  }

} // end namespace Antioch

#endif // ANTIOCH_DATA_EQUILIBRIUM_H
