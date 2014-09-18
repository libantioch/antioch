//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_MOLECULAR_BINARY_DIFFUSION_UTILS_H
#define ANTIOCH_MOLECULAR_BINARY_DIFFUSION_UTILS_H

#include "antioch/molecular_binary_diffusion_utils_decl.h"

namespace Antioch
{
   // getting tag
   template <typename CoeffType, typename Interpolator>
   struct physical_tag_type<MolecularBinaryDiffusion<CoeffType,Interpolator> >
   {
      typedef bimolecular_diffusion_tag type;
   };

   // physical set boolean
   template<typename CoeffType, typename Interpolator>
   struct is_physical_set<MolecularBinaryDiffusion<CoeffType, Interpolator> >
   {
      static const bool value = true;
   };

   template <typename CoeffType, typename Interpolator>
   struct SetOrEquation<MolecularBinaryDiffusion<CoeffType,Interpolator>,true>
   {
      typedef std::vector<std::vector<BinaryDiffusion<CoeffType,Interpolation> >* > type;
   }

   // we can initialize without the user's help,
   // so we tag the physical_set_tag:
   // we need to define 
   // physical_set_initialize
   // physical_set_delete
   // physical_set_add
   // physical_set_reset
   // but not the default ones that are not concerned
   template <typename CoeffType, typename Interpolator>
   struct physical_set_tag<MolecularBinaryDiffusion<CoeffType,Interpolator> >
   {
     typedef typename molecular_binary_diffusion_tag type;
   };

   // we can initialize without the user's help
   template <typename ModelSet>
   void physical_set_initialize(ModelSet & mod, molecular_binary_diffusion_tag )
   {
      build_molecular_binary_diffusion(mod);
   }

   template <typename ModelSet>
   void physical_set_delete(ModelSet & mod, molecular_binary_diffusion_tag )
   {}

   template <typename Model, typename InitType>
   void physical_set_add(unsigned int s, SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, molecular_binary_diffusion_tag)
   {
     antioch_assert(!set[s]);
     set[s] = new Model(init);
   }

   template <typename Model, typename InitType>
   void physical_set_rest(unsigned int s, SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, molecular_binary_diffusion_tag)
   {
     set[s]->reset_coeffs(init);
   }



   template<typename Model, typename StateType, typename MatrixStateType>
   void physical_set_operator_diffusion_comes_first(const Model & set, const StateType & T, const StateType & cTot, MatrixStateType & Ds, bimolecular_diffusion_tag)
   {
       antioch_assert_equal_to(Ds.size(),set.size());

       for(unsigned int i = 0; i < Ds.size(); i++)
       {
         antioch_assert_equal_to(Ds[i].size(),set.size());
         for(unsigned int j = 0; j <= i; j++)
         {
             Ds[i][j] = _set[i][j](T,cTot);
             Ds[j][i] = Ds[i][j];
         }
       }
   }

   template<typename Model, typename StateType, typename VectorStateType>
   void physical_set_operator_diffusion_comes_first(const Model & set, const StateType & T, const StateType & cTot, VectorStateType & Ds, bimolecular_diffusion_second_tag)
   {
       antioch_assert_equal_to(Ds.size(),set.size());

       for(unsigned int i = 0; i < Ds.size(); i++)
       {
             Ds[i] = _set[i][i](T,cTot);
       }
   }

        // template around mixture mainly to avoid 
        // ChemicalMixture<CoeffType> forward declaration
        // Formulae used is
        // d_s = 1 - y_s / (sum_{j \neq s} x_j/D_{sj})
   template <typename Mixture, typename MatrixStateType, typename VectorStateType>
   void wilke_diffusion_rule(const Mixture & mixture, const VectorStateType & mass_fractions, const MatrixStateType & Ds, VectorStateType & ds,bimolecular_diffusion_tag)
   {
       antioch_assert(ds.size(),mixture.n_species());
       antioch_assert(ds.size(),mass_fractions.size());

       VectorStateType molar_fractions = zero_clone(mass_fractions);
       _transport_mixture.chemical_mixture().X(_transport_mixture.chemical_mixture().M(mass_fractions),mass_fractions,molar_fractions);

       for(unsigned int s = 0; s < ds.size(); s++)
       {
          ds[s] = constant_clone(mass_fractions[s],1) - mass_fractions[s];
          typename value_type<VectorStateType>::type denom = zero_clone(molar_mass[0]);
          for(unsigned int j = 0; j < ds.size(); j++)
          {
             if(j == s)continue;
             denom += molar_fractions[j] / Ds[s][j];
          }
          ds[s] /= denom;
       }
   }

}

#endif
