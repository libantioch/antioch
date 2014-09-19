//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_MOLECULAR_BINARY_DIFFUSION_UTILS_H
#define ANTIOCH_MOLECULAR_BINARY_DIFFUSION_UTILS_H

#include "antioch/transport_species.h"
#include "antioch/molecular_binary_diffusion_utils_decl.h"

namespace Antioch
{
   // getting tag
   template <typename CoeffType, typename Interpolator>
   struct physical_tag<MolecularBinaryDiffusion<CoeffType,Interpolator> >
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
     typedef default_physical_tag      viscosity_type;
     typedef bimolecular_diffusion_tag diffusion_species_type;
     typedef default_physical_tag      diffusion_mixture_type;
     typedef default_physical_tag      thermal_conductivity_type;
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
      typedef std::vector<std::vector<BinaryDiffusion<CoeffType,Interpolation>* > > type;
   }

   // initializer
  template<typename CoeffType, typename Interpolator>
  struct Initializer<MolecularBinaryDiffusion<CoeffType,Interpolator> >
  {
     unsigned int j;
     TransportSpecies<CoeffType> & si;
     TransportSpecies<CoeffType> & sj;
  };

   // we can initialize without the user's help
   template <typename ModelSet>
   void physical_set_initialize(ModelSet & mod, molecular_binary_diffusion_tag )
   {
      mod.set().resize(mod.mixture().n_species());
      for(unsigned int i = 0; i < mod.set().size(); i++)
      {
          mod.set()[i].resize(i,NULL);
      }
      build_molecular_binary_diffusion(mod);
   }


   template <typename Model, typename InitType>
   void physical_set_add(Model & set, const InitType & init, bimolecular_diffusion_tag){}

   template <typename Model, typename InitType>
   void physical_set_add(unsigned int s, typename SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, molecular_binary_diffusion_tag)
   {
     antioch_assert_less(s,set.size());
     antioch_assert_less(init.j,set[s].size());
     antioch_assert(!set[s][init.j]);

     set[s][init.j] = new Model(init.si,init.sj);
   }


   template <typename Model, typename InitType>
   void physical_set_reset(Model & set, const InitType & init, bimolecular_diffusion_tag)
   {}

   template <typename Model, typename InitType>
   void physical_set_reset(unsigned int s, SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const InitType & init, molecular_binary_diffusion_tag)
   {
     set[s][init.j]->reset_coeffs(init.si,init.sj);
   }

   template <typename Model>
   void physical_set_print(SetOrEquation<Model,is_physical_set<Model>::value>::type & set, const std::map<unsigned int, std::string>& spec, std::ostream & out, bimolecular_diffusion_tag)
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


   template<typename Model, typename StateType, typename MatrixStateType>
   void physical_set_operator_diffusion(const Model & set, const StateType & T, const StateType & cTot, MatrixStateType & Ds, bimolecular_diffusion_tag)
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

   template<typename Model, typename StateType>
   void physical_set_operator_diffusion(unsigned int s, const Model & set, const StateType & T, const StateType & cTot, StateType & Ds, bimolecular_diffusion_second_tag)
   {
       antioch_assert_equal_to(Ds.size(),set.size());
       Ds = _set[s][s](T,cTot);
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
       _mixture.X(_mixture.M(mass_fractions),mass_fractions,molar_fractions);

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
