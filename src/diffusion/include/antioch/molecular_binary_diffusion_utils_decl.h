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

   // getting tag
   template <typename CoeffType, typename Interpolator>
   struct physical_tag<MolecularBinaryDiffusion<CoeffType,Interpolator> >;

   // physical set boolean
   template<typename CoeffType, typename Interpolator>
   struct is_physical_set<MolecularBinaryDiffusion<CoeffType, Interpolator> >;

   // custom set
   template <typename CoeffType, typename Interpolator, bool B>
   struct SetOrEquation<MolecularBinaryDiffusion<CoeffType,Interpolator>,B>;

   // custom initialization
   template<typename CoeffType, typename Interpolator>
   struct Initializer<MolecularBinaryDiffusion<CoeffType,Interpolator> >;

   // we can initialize without the user's help
   template <typename ModelSet>
   void physical_set_initialize(ModelSet & mod, bimolecular_diffusion_tag );

   // custom delete
   template <typename ModelSet>
   void physical_set_delete(ModelSet & mod, bimolecular_diffusion_tag );

   // custom add
   template <typename Model, typename InitType>
   void physical_set_add(unsigned int s, typename SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, bimolecular_diffusion_tag);

   // requires to void this
   template <typename Model, typename InitType>
   void physical_set_add(Model & set, const InitType & init, bimolecular_diffusion_tag);

   // custom reset
   template <typename Model, typename InitType>
   void physical_set_reset(unsigned int s, typename SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, bimolecular_diffusion_tag);

   // requires to void this
   template <typename Model, typename InitType>
   void physical_set_reset(Model & set, const InitType & init, bimolecular_diffusion_tag);

   // custom print
   template <typename Model>
   void physical_set_print(SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const std::map<unsigned int, std::string>& spec, std::ostream & out, bimolecular_diffusion_tag);

   template<typename Model, typename StateType, typename MatrixStateType>
   void physical_set_operator_diffusion(const Model & set, const StateType & T, const StateType & cTot, MatrixStateType & Ds, bimolecular_diffusion_tag);

   template<typename Model, typename StateType>
   void physical_set_operator_diffusion(unsigned int s, const Model & set, const StateType & T, const StateType & cTot, StateType & Ds, bimolecular_diffusion_tag);

        // template around mixture mainly to avoid 
        // ChemicalMixture<CoeffType> forward declaration
        // Formulae used is
        // d_s = 1 - y_s / (sum_{j \neq s} x_j/D_{sj})
   template <typename Mixture, typename MatrixStateType, typename VectorStateType>
   void wilke_diffusion_rule(const Mixture & mixture, const VectorStateType & mass_fractions, const MatrixStateType & Ds, VectorStateType & ds, bimolecular_diffusion_tag);
}

#endif
