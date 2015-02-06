//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014 Paul T. Bauman, Benjamin S. Kirk, Sylvain Plessis,
//                    Roy H. Stonger
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

#ifndef ANTIOCH_MOLECULAR_BINARY_DIFFUSION_UTILS_H
#define ANTIOCH_MOLECULAR_BINARY_DIFFUSION_UTILS_H

#include "antioch/transport_species.h"
#include "antioch/molecular_binary_diffusion_utils_decl.h"
#include "antioch/molecular_binary_diffusion_building.h"

namespace Antioch
{
   // getting tag
   template <typename CoeffType, typename Interpolator>
   struct physical_tag<MolecularBinaryDiffusion<CoeffType,Interpolator> >:
        public physical_tag_base<MolecularBinaryDiffusion<CoeffType,Interpolator> >
   {
      typedef bimolecular_diffusion_tag type;
        // kind of set tag
     typedef bimolecular_diffusion_tag set_type;
        // some models require specific initialization
        // (typically automatic initialization)
     typedef bimolecular_diffusion_tag init_type;
        // some models require specific initialization
        // but not specific deletion
     typedef bimolecular_diffusion_tag del_type;
        // for operators, diffusion is special, see comment below
     typedef bimolecular_diffusion_tag diffusion_species_type;
   };

   // physical set boolean
   template<typename CoeffType, typename Interpolator>
   struct is_physical_set<MolecularBinaryDiffusion<CoeffType, Interpolator> >
   {
      static const bool value = true;
   };

        // special set
   template <typename CoeffType, typename Interpolator>
   struct SetOrEquation<MolecularBinaryDiffusion<CoeffType,Interpolator>,true>
   {
      typedef std::vector<std::vector<MolecularBinaryDiffusion<CoeffType,Interpolator>* > > type;
   };

   // initializer
  template<typename CoeffType, typename Interpolator>
  struct Initializer<MolecularBinaryDiffusion<CoeffType,Interpolator>, bimolecular_diffusion_tag >
  {
      Initializer(unsigned int i, const TransportSpecies<CoeffType> & s1, const TransportSpecies<CoeffType> & s2):
                j(i),si(s1),sj(s2){}

     unsigned int j;
     const TransportSpecies<CoeffType> & si;
     const TransportSpecies<CoeffType> & sj;
  };

   // we can initialize without the user's help
   template <typename Model>
   void physical_set_initialize(Model & mod, bimolecular_diffusion_tag )
   {
      mod.set().resize(mod.mixture().n_species());
      for(unsigned int i = 0; i < mod.set().size(); i++)
      {
          mod.set()[i].resize(i+1,NULL);
      }
      build_molecular_binary_diffusion(mod);
   }

   // we can initialize without the user's help
   template <typename Model>
   void physical_set_delete(Model & mod, bimolecular_diffusion_tag )
   {
      for(unsigned int i = 0; i < mod.set().size(); i++)
      {
        for(unsigned int j = 0; j < mod.set()[i].size(); j++)
        {
          if(mod.set()[i][j])delete mod.set()[i][j];
          mod.set()[i][j] = NULL;
        }
      }
   }



   template <typename Model, typename InitType>
   void physical_set_add(Model & set, const InitType & init, bimolecular_diffusion_tag){}

   template <typename Model, typename InitType>
   void physical_set_add(unsigned int s, typename SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, bimolecular_diffusion_tag)
   {
     antioch_assert_less(s,set.size());
     antioch_assert_less(init.j,set[s].size());
     antioch_assert(!set[s][init.j]);

        // matrix is symmetric
     (s >= init.j)?set[s][init.j] = new Model(init.si,init.sj):
                   set[init.j][s] = new Model(init.sj,init.si);
   }


   template <typename Model, typename InitType>
   void physical_set_reset(Model & set, const InitType & init, bimolecular_diffusion_tag)
   {}

   template <typename Model, typename InitType>
   void physical_set_reset(unsigned int s, typename SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, bimolecular_diffusion_tag)
   {
        // matrix is symmetric
     (s >= init.j)?set[s][init.j]->reset_coeffs(init.si,init.sj):
                   set[init.j][s]->reset_coeffs(init.sj,init.si);
   }

   template <typename CoeffType, typename Interpolator>
   void physical_set_print(typename SetOrEquation<MolecularBinaryDiffusion<CoeffType,Interpolator>,true>::type & set, 
                            const std::map<unsigned int, std::string>& spec, std::ostream & out, bimolecular_diffusion_tag)
   {
        out << "Bimolecular diffusion" << std::endl;

        for(unsigned int s = 0; s < spec.size(); ++s)
        {
            for(unsigned int j = 0; j <= s; ++j)
            {
                out << spec.at(s) << ":" << spec.at(j) << " " << *set[s][j] << std::endl;
            }
        }
   }


      // here we need only the lower triangular
      // matrix, the diagonal need not be computed
      // for the diffusion, it is required only
      // for the thermal conduction in case of
      // kinetics theory
   template<typename Model, typename StateType, typename MatrixStateType>
   void physical_set_operator_diffusion(const Model & set, const StateType & T, const StateType & cTot, MatrixStateType & Ds, bimolecular_diffusion_tag)
   {
       antioch_assert_equal_to(Ds.size(),set.size());

       for(unsigned int i = 0; i < Ds.size(); i++)
       {
         antioch_assert_equal_to(Ds[i].size(),set.size());
         for(unsigned int j = 0; j < i; j++)
         {
             Ds[i][j] = (*set[i][j])(T,cTot);
         }
       }
   }

     // self-diffusion
   template<typename Model, typename StateType>
   void physical_set_operator_diffusion(unsigned int s, const Model & set, const StateType & T, const StateType & cTot, StateType & Ds, bimolecular_diffusion_tag)
   {
       Ds = (*set[s][s])(T,cTot);
   }

        // template around mixture mainly to avoid 
        // ChemicalMixture<CoeffType> forward declaration
        // Formulae used is
        // d_s = (1 - y_s) / (sum_{j \neq s} x_j/D_{sj})
        // might want to change it for stability
   template <typename Mixture, typename MatrixStateType, typename VectorStateType>
   void wilke_diffusion_rule(const Mixture & mixture, const VectorStateType & mass_fractions, const MatrixStateType & Ds, VectorStateType & ds,bimolecular_diffusion_tag)
   {
       antioch_assert_equal_to(ds.size(),mixture.n_species());
       antioch_assert_equal_to(ds.size(),mass_fractions.size());

       VectorStateType molar_fractions = zero_clone(mass_fractions);
       mixture.X(mixture.M(mass_fractions),mass_fractions,molar_fractions);

       for(unsigned int s = 0; s < ds.size(); s++)
       {
          ds[s] = constant_clone(mass_fractions[s],1) - mass_fractions[s];
          typename value_type<VectorStateType>::type denom = zero_clone(mass_fractions[0]);
          for(unsigned int j = 0; j < ds.size(); j++)
          {
             if(j == s)continue;
             (s > j)?denom += molar_fractions[j] / Ds[s][j]:
                     denom += molar_fractions[j] / Ds[j][s];
          }
          ds[s] /= denom;
       }
   }

}

#endif
