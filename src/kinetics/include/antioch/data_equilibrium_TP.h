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

#ifndef ANTIOCH_DATA_EQUILIBRIUM_TP_H
#define ANTIOCH_DATA_EQUILIBRIUM_TP_H

// antioch
#include "antioch/data_equilibrium.h"
#include "antioch/chemical_mixture.h"

//C++
#include <vector>

namespace Antioch
{

//forward declaration
  template<typename CoeffType>
  class ChemicalMixture;

  template<typename CoeffType = double>
  class DataEquilibriumTP: public DataEquilibrium<CoeffType>
  {
     public:
          DataEquilibriumTP(const CoeffType &T_mix, const CoeffType &P_mix, 
                    const ReactionSet<CoeffType> &reac_set);
          ~DataEquilibriumTP();

          void fill_constrain(const std::vector<CoeffType> &molar_densities,
                              std::vector<CoeffType> &F,std::vector<std::vector<CoeffType> > &jacob);//constant pressure


          const CoeffType local_pressure(const std::vector<CoeffType> &molar_densities) const {return loc_P;}

          const unsigned int n_constrain() const {return 1;}
          void update_constrain(unsigned int icstr, const CoeffType & delta_cstr);//constrain on P

     private:
     CoeffType loc_P;

//I don't want it to be used, change it if you want
     DataEquilibriumTP(){}
  };


  template<typename CoeffType>
  inline
  DataEquilibriumTP<CoeffType>::DataEquilibriumTP(const CoeffType &T_mix, const CoeffType &P_mix, 
                    const ReactionSet<CoeffType> &reac_set):
     DataEquilibrium<CoeffType>(T_mix,P_mix,reac_set)
  {
     return;
  }

  template<typename CoeffType>
  inline
  DataEquilibriumTP<CoeffType>::~DataEquilibriumTP()
  {
    return;
  }

  template<typename CoeffType>
  inline
  void DataEquilibriumTP<CoeffType>::fill_constrain(const std::vector<CoeffType> &molar_densities,std::vector<CoeffType> &F, std::vector<std::vector<CoeffType> > &jacob)
  {

     antioch_assert_equal_to(F.size(),this->_reac_set.n_species());
     antioch_assert_equal_to(jacob.size(),this->_reac_set.n_species());
     for(unsigned int i = 0; i < this->_reac_set.n_species(); i++)
     {
       antioch_assert_equal_to(jacob[i].size(),this->_reac_set.n_species());
     }

     CoeffType total_dens = Antioch::zero_clone(this->_P);
     CoeffType molar_sum  = Antioch::zero_clone(this->_P);
     total_dens = loc_P/(Constants::R_universal<CoeffType>()*this->_T);
     jacob.push_back(std::vector<CoeffType>());
     jacob[this->_reac_set.n_species()].resize(this->_reac_set.n_species()+1);
     for(unsigned int i = 0; i < this->_reac_set.n_species(); i++)
     {
        molar_sum += molar_densities[i];
        jacob[this->_reac_set.n_species()][i] = 1./this->_reac_set.chemical_mixture().M(i);
        jacob[i].push_back(0.);
     }
     jacob[this->_reac_set.n_species()][this->_reac_set.n_species()] = -1.;
     F.push_back(molar_sum - total_dens);
  }

  template<typename CoeffType>
  inline
  void DataEquilibriumTP<CoeffType>::update_constrain(unsigned int icstr, const CoeffType &delta_cstr)
  {
     antioch_assert_equal_to(icstr,0);
     loc_P += delta_cstr * Constants::R_universal<CoeffType>() * this->_T;
  }

} // end namespace Antioch

#endif // ANTIOCH_DATA_EQUILIBRIUM_H
