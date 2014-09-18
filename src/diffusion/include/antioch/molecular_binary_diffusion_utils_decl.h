//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_MOLECULAR_BINARY_DIFFUSION_UTILS_DECL_H
#define ANTIOCH_MOLECULAR_BINARY_DIFFUSION_UTILS_DECL_H

#include "antioch/physics_metaprogramming_decl.h"

namespace Antioch
{
// binary molecular
   // tag
   struct bimolecular_diffusion_tag{};
   struct bimolecular_diffusion_second_tag{};

   // getting tag
   template <typename CoeffType, typename Interpolator>
   struct physical_tag_type<MolecularBinaryDiffusion<CoeffType,Interpolator> >;

   // physical set boolean
   template<typename CoeffType, typename Interpolator>
   struct is_physical_set<MolecularBinaryDiffusion<CoeffType, Interpolator> >;

   template <typename CoeffType, typename Interpolator, bool B>
   struct SetOrEquation<MolecularBinaryDiffusion<CoeffType,Interpolator>,B>;

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
     typedef typename pure_species_viscosity_tag type;
   };

   // we can initialize without the user's help
   template <typename ModelSet>
   void physical_set_initialize(ModelSet & mod, pure_species_viscosity_tag )
   {
      build_pure_species_viscosity(mod);
   }

   template <typename ModelSet>
   void physical_set_delete(ModelSet & mod, pure_species_viscosity_tag )
   {}

   template <typename Model, typename InitType>
   void physical_set_add(unsigned int s, SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, pure_species_viscosity_tag)
   {
     antioch_assert(!set[s]);
     set[s] = new Model(init);
   }

   template <typename Model, typename InitType>
   void physical_set_rest(unsigned int s, SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, pure_species_viscosity_tag)
   {
     set[s]->reset_coeffs(init);
   }




   template<typename Model, typename StateType, typename MatrixStateType>
   void physical_set_operator_diffusion_comes_first(const Model & set, const StateType & T, const StateType & cTot, MatrixStateType & Ds, bimolecular_diffusion_tag);

   template<typename Model, typename StateType, typename VectorStateType>
   void physical_set_operator_diffusion_comes_first(const Model & set, const StateType & T, const StateType & cTot, VectorStateType & Ds, bimolecular_diffusion_second_tag);

        // template around mixture mainly to avoid 
        // ChemicalMixture<CoeffType> forward declaration
        // Formulae used is
        // d_s = 1 - y_s / (sum_{j \neq s} x_j/D_{sj})
   template <typename Mixture, typename MatrixStateType, typename VectorStateType>
   void wilke_diffusion_rule(const Mixture & mixture, const VectorStateType & mass_fractions, const MatrixStateType & Ds, VectorStateType & ds,bimolecular_diffusion_tag);
}

#endif
