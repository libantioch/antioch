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
#include <string>

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
                    const ReactionSet<CoeffType> &reac_set, const std::string &key = std::string());
          virtual ~DataEquilibrium();

          void fill_constrain(const std::vector<CoeffType> &molar_densities,
                                        std::vector<CoeffType> &F,std::vector<std::vector<CoeffType> > &jacob);

          void set_constrain(const std::string &key);

          const CoeffType T() const {return _T;}
          const CoeffType P() const {return _P;}
          const CoeffType local_pressure(const std::vector<CoeffType> &molar_densities) const; 
          const ReactionSet<CoeffType> &reaction_set() const {return _reac_set;}

          const unsigned int n_constrain() const {return _n_constrain;}
          void update_constrain(unsigned int icstr, const CoeffType &delta_cstr);

          void set_mass(const CoeffType &m) {_fixed_mass = m;}
          void set_T(const CoeffType &t)    {_T = t;}
          void set_P(const CoeffType &p)    {_P = p;}


     protected:
     const CoeffType &_T;
     const CoeffType &_P;
     const ReactionSet<CoeffType> &_reac_set;

     private:
//I don't want it to be used, change it if you want
     DataEquilibrium(){return;}


     CoeffType _fixed_pressure;
     CoeffType _fixed_mass;
     unsigned int m_constrain;
     unsigned int a_constrain;
     unsigned int p_constrain;
     unsigned int _n_constrain;

  };


  template<typename CoeffType>
  inline
  DataEquilibrium<CoeffType>::DataEquilibrium(const CoeffType &T_mix, const CoeffType &P_mix, 
                    const ReactionSet<CoeffType> &reac_set, const std::string &key):
     _T(T_mix),_P(P_mix),_reac_set(reac_set),m_constrain(1),a_constrain(0),p_constrain(0),_n_constrain(0)
  {
     _fixed_pressure = this->_P/(Constants::R_universal<CoeffType>() * this->_T);
     set_constrain(key);
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
  void DataEquilibrium<CoeffType>::set_constrain(const std::string &key)
  {
     if(key.find("m") != std::string::npos)m_constrain = 1;
     if(key.find("a") != std::string::npos)a_constrain = 1;
     if(key.find("p") != std::string::npos)p_constrain = 1;
     _n_constrain = key.length();
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

  template<typename CoeffType>
  inline
  void DataEquilibrium<CoeffType>::fill_constrain(const std::vector<CoeffType> &molar_densities,std::vector<CoeffType> &F, std::vector<std::vector<CoeffType> > &jacob)
  {
     if(m_constrain)
     {
       antioch_assert_equal_to(F.size(),jacob.size());
       for(unsigned int i = 0; i < F.size(); i++)
       {
         antioch_assert_equal_to(jacob[i].size(),F.size());
       }

       CoeffType mass_sum  = Antioch::zero_clone(this->_P);
       for(unsigned int i = 0; i < this->_reac_set.n_species()-1; i++)
       {
          mass_sum += molar_densities[i] * this->_reac_set.chemical_mixture().M(i);
          jacob.back()[i] = 1.;
       }
       jacob.back().back() = 1.;
       F.back() = mass_sum - _fixed_mass;
     }
     if(p_constrain)
     {
       antioch_assert_equal_to(F.size(),jacob.size());
       for(unsigned int i = 0; i < F.size(); i++)
       {
         antioch_assert_equal_to(jacob[i].size(),F.size());
       }

       CoeffType molar_sum  = Antioch::zero_clone(this->_P);
       jacob.push_back(std::vector<CoeffType>(F.size()+1,0.));
       for(unsigned int i = 0; i < this->_reac_set.n_species(); i++)
       {
          molar_sum += molar_densities[i];
          jacob.back()[i] = 1./this->_reac_set.chemical_mixture().M(i);
          jacob[i].push_back(0.);
       }
       jacob.back().back() = -1.;
       F.push_back(molar_sum - _fixed_pressure);
     }
  }

} // end namespace Antioch

#endif // ANTIOCH_DATA_EQUILIBRIUM_H
