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
#ifndef ANTIOCH_EQUILIBRIUM_EVALUATOR_H
#define ANTIOCH_EQUILIBRIUM_EVALUATOR_H

#include "antioch_config.h"

#ifdef ANTIOCH_HAVE_EIGEN
#include "antioch/eigen_utils_decl.h"
#include "antioch/eigen_utils.h"

// antioch
#include "antioch/kinetics_evaluator.h"
#include "antioch/cea_thermo.h"
#include "antioch/data_equilibrium.h"


//C++
#include <vector>

namespace Antioch
{

  template<typename CoeffType = double,typename StateType=CoeffType>
  class EquilibriumEvaluator
  {
  public:
    EquilibriumEvaluator(DataEquilibrium<CoeffType> &data_eq,
                         KineticsEvaluator<CoeffType,StateType>&kin,
                         const long double tol = std::numeric_limits<StateType>::epsilon() * 100);

    ~EquilibriumEvaluator();

    void first_guess_mass_fraction(const std::vector<CoeffType> &first) {eq_mass_fraction = first;}
    void first_guess_molar_fraction(const std::vector<CoeffType> &first);


    const std::vector<CoeffType> mass_fraction_equilibrium()   const {return eq_mass_fraction;}
    const std::vector<CoeffType> molar_densities_equilibrium() const {return eq_molar_densities;}
    const CoeffType Peq() const;

    void equilibrium();

    const bool success()              const {return !over_threshold;}
    const double residual()           const {return thres;}
    const double conv_threshold()     const {return threshold;}
    const unsigned int max_loop_tol() const {return max_loop;}

    void provide_first_guess(const std::vector<CoeffType> &first_mass_fraction)
    {eq_mass_fraction = first_mass_fraction;}


  private:

    void first_guess_iso();
    void init_from_mass();
    void calculate_function_and_jacobian(std::vector<StateType> &F,
                                         std::vector<std::vector<StateType> > &jacob);

    std::vector<CoeffType> eq_molar_densities;
    std::vector<CoeffType> eq_mass_fraction;
    std::vector<CoeffType> eq_mass;
    CoeffType mass_tot;
    CoeffType mass_tot_ini;
    KineticsEvaluator<CoeffType,StateType> &kinetics_eval;
    DataEquilibrium<CoeffType> & data_storage_and_constrain;


    //useless but needed
    std::vector<CoeffType> dGibbs_RT_dT;
    std::vector<CoeffType> dmass_dT;

    //local variables
    const long double threshold;
    const unsigned int max_loop;
    bool over_threshold;
    double thres;

    //! 
    EquilibriumEvaluator();
  };

  template<typename CoeffType, typename StateType>
  inline
  EquilibriumEvaluator<CoeffType,StateType>::~EquilibriumEvaluator()
  {return;}

  template<typename CoeffType, typename StateType>
  inline
  EquilibriumEvaluator<CoeffType,StateType>::EquilibriumEvaluator( DataEquilibrium<CoeffType>& data_eq,
                                                                   KineticsEvaluator<CoeffType,StateType>& kin,
                                                                   const long double tol )
    : mass_tot(0.),
      kinetics_eval(kin),
      data_storage_and_constrain(data_eq),
      threshold(tol),max_loop(100),
      over_threshold(true),thres(-1.)
  {
    eq_mass_fraction.resize(data_storage_and_constrain.reaction_set().n_species(),0.);
    first_guess_iso();
    return;
  }

  template<typename CoeffType, typename StateType>
  inline
  void EquilibriumEvaluator<CoeffType, StateType>::equilibrium()
  {

    init_from_mass();

    std::vector<StateType> F;
    std::vector<std::vector<StateType> > jacob;
    unsigned int ncol(data_storage_and_constrain.reaction_set().n_species());
    unsigned int nrow(data_storage_and_constrain.reaction_set().n_species());

    Eigen::Matrix<StateType,Eigen::Dynamic,Eigen::Dynamic> jacob_decomp(nrow,ncol);
    Eigen::Matrix<StateType,Eigen::Dynamic,1> deltaX(nrow);
    Eigen::Matrix<StateType,Eigen::Dynamic,1> minus_function(nrow);

    over_threshold = true;
    unsigned int nloop(0);
    while(over_threshold)
      {
        nloop++;
        calculate_function_and_jacobian(F,jacob);
        Antioch::set_zero(thres);
        antioch_assert_equal_to(F.size(),ncol);
        antioch_assert_equal_to(jacob.size(),ncol);
        for(unsigned int i = 0; i < nrow; i++)
          {
            thres += (F[i] > 0.)?F[i]:-F[i];
            minus_function(i) = -F[i];
            antioch_assert_equal_to(jacob[i].size(),F.size());
            for(unsigned int j = 0; j < ncol; j++)
              {
                jacob_decomp(i,j) = jacob[i][j];
              }
          }
        over_threshold = thres > threshold;
        if(!over_threshold)break;

        std::cout << jacob_decomp << " and\n " << minus_function << std::endl;
        deltaX = jacob_decomp.partialPivLu().solve(minus_function);
        std::cout << "gives\n" << deltaX << std::endl;

        Antioch::set_zero(thres);
        Antioch::set_zero(mass_tot);
        for(unsigned int i = 0; i < data_storage_and_constrain.reaction_set().n_species(); i++)
          {
            thres += (deltaX(i) > 0.)?deltaX(i):-deltaX(i);

            //std::cout << data_storage_and_constrain.reaction_set().chemical_mixture().chemical_species()[i]->species() << " (" << minus_function(i) << "): " 
            //          << eq_mass[i] << " ** ";

            eq_mass[i] += deltaX(i);
            if(eq_mass[i] < 0.)eq_mass[i] = -eq_mass[i];
            eq_molar_densities[i] = eq_mass[i]/data_storage_and_constrain.reaction_set().chemical_mixture().M(i);
            //std::cout << eq_mass[i] <<  std::endl;
            /*
              eq_molar_densities[i] += deltaX(i);
              eq_mass[i] = eq_molar_densities[i] * data_storage_and_constrain.reaction_set().chemical_mixture().M(i);
              std::cout << eq_molar_densities[i] <<  std::endl;
            */
            mass_tot += eq_mass[i];
          }

        if(thres != thres)break;
        over_threshold = thres > threshold;
        if(nloop > max_loop)break;
      }

  }


  template<typename CoeffType, typename StateType>
  inline
  void EquilibriumEvaluator<CoeffType, StateType>::calculate_function_and_jacobian(std::vector<StateType>& F, std::vector<std::vector<StateType> > &jacob)
  {
    F = Antioch::zero_clone(eq_mass_fraction);//first fill with pure kinetics
    ///kinetics part only first
    //updating
    jacob.resize(data_storage_and_constrain.reaction_set().n_species());

    for(unsigned int i = 0; i < data_storage_and_constrain.reaction_set().n_species(); i++)
      {
        jacob[i].resize(data_storage_and_constrain.reaction_set().n_species(),0.);
        eq_mass_fraction[i] = eq_mass[i]/mass_tot;
      }

    //setting the system
    std::vector<CoeffType> h_RT_minus_s_R(data_storage_and_constrain.reaction_set().n_species());

    Antioch::CEAThermodynamics<StateType> thermo(data_storage_and_constrain.reaction_set().chemical_mixture() );
    typedef typename Antioch::CEAThermodynamics<StateType>::template Cache<StateType> Cache;
    thermo.h_RT_minus_s_R(Cache(data_storage_and_constrain.T()),h_RT_minus_s_R);

    //omega_dot & jacobian calculations
    /*std::cout << "T: " << data_storage_and_constrain.T()<< std::endl;
      std::cout << "local P: " << cur_P << std::endl;
      std::cout << "rho: " << rho << std::endl;
      std::cout << "R_mix: " << R_mix << std::endl;
      std::cout << "n_species: " << data_storage_and_constrain.reaction_set().n_species() << std::endl;
      std::cout << "n_reaction: " << data_storage_and_constrain.reaction_set().n_reactions() << std::endl;
    */   kinetics_eval.compute_mass_sources_and_derivs(data_storage_and_constrain.T(),
                                                       eq_molar_densities,
                                                       h_RT_minus_s_R,
                                                       dGibbs_RT_dT,
                                                       F,
                                                       dmass_dT,
                                                       jacob);

    /*  std::cout << "F size " << F.size() << std::endl << "\t";
        for(unsigned int i = 0; i < F.size(); i++)std::cout << F[i] << " ";
        std::cout << std::endl;
    */
    /*
      for(unsigned int i = 0; i < jacob.size(); i++)
      {
      F[i] /= data_storage_and_constrain.reaction_set().chemical_mixture().M(i);
      for(unsigned int j= 0; j < jacob[i].size(); j++)
      {
      jacob[i][j] /= data_storage_and_constrain.reaction_set().chemical_mixture().M(i)/data_storage_and_constrain.reaction_set().chemical_mixture().M(j);
      }
      }
    */

    data_storage_and_constrain.fill_constrain(eq_molar_densities,F,jacob); //constrain added
  }

  template<typename CoeffType, typename StateType>
  inline
  void EquilibriumEvaluator<CoeffType, StateType>::first_guess_iso()
  {
    CoeffType iso = 1./(double)(data_storage_and_constrain.reaction_set().n_species());
    std::fill(eq_mass_fraction.begin(),eq_mass_fraction.end(),iso);
  }


  template<typename CoeffType, typename StateType>
  inline
  void EquilibriumEvaluator<CoeffType, StateType>::init_from_mass()
  {
    eq_molar_densities = Antioch::zero_clone(eq_mass_fraction);
    eq_mass = Antioch::zero_clone(eq_mass_fraction);
    dGibbs_RT_dT = Antioch::zero_clone(eq_mass_fraction);
    dmass_dT = Antioch::zero_clone(eq_mass_fraction);
    mass_tot_ini = data_storage_and_constrain.reaction_set().chemical_mixture().M(eq_mass_fraction) *
      data_storage_and_constrain.P() /
      (Constants::R_universal<CoeffType>() * data_storage_and_constrain.T());
    mass_tot = mass_tot_ini;
    data_storage_and_constrain.set_mass(mass_tot);

    for(unsigned int i = 0; i < eq_mass_fraction.size(); i++)
      {
        eq_mass[i] = mass_tot * eq_mass_fraction[i];
        eq_molar_densities[i] = eq_mass[i] / data_storage_and_constrain.reaction_set().chemical_mixture().M(i);
      }
  }

  template<typename CoeffType, typename StateType>
  inline
  const CoeffType EquilibriumEvaluator<CoeffType, StateType>::Peq() const
  {
    return mass_tot / (data_storage_and_constrain.reaction_set().chemical_mixture().M(eq_mass_fraction)) *
      (Constants::R_universal<CoeffType>()*data_storage_and_constrain.T());
  }



  template<typename CoeffType, typename StateType>
  inline
  void EquilibriumEvaluator<CoeffType, StateType>::first_guess_molar_fraction(const std::vector<CoeffType> &first)
  {
    antioch_assert_equal_to(first.size(),data_storage_and_constrain.reaction_set().n_species() );
    CoeffType denom;
    Antioch::set_zero(denom);
    for(unsigned int i = 0; i < data_storage_and_constrain.reaction_set().n_species(); i++)
      {
        denom += first[i] * data_storage_and_constrain.reaction_set().chemical_mixture().M(i);
      }
    for(unsigned int i = 0; i < data_storage_and_constrain.reaction_set().n_species(); i++)
      {
        eq_mass_fraction[i] = first[i] * data_storage_and_constrain.reaction_set().chemical_mixture().M(i)/denom;
      }

    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_EQUILIBRIUM_EVALUATOR_H
#endif //ANTIOCH_HAVE_EIGEN
