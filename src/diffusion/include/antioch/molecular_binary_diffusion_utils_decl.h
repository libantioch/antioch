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
