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
   struct Initializer<MolecularBinaryDiffusion<CoeffType,Interpolator>, bimolecular_diffusion_tag >;

   // we can initialize without the user's help
   template <typename Model>
   void physical_set_initialize(Model & mod, bimolecular_diffusion_tag );

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
   template <typename CoeffType, typename Interpolator>
   void physical_set_print(typename SetOrEquation<MolecularBinaryDiffusion<CoeffType,Interpolator>,true>::type & set, 
                            const std::map<unsigned int, std::string>& spec, std::ostream & out, bimolecular_diffusion_tag);

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
