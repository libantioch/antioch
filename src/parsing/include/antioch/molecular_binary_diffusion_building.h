//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_MOLECULAR_BINARY_DIFFUSION_BUILDING_H
#define ANTIOCH_MOLECULAR_BINARY_DIFFUSION_BUILDING_H

// Antioch
#include "antioch/molecular_binary_diffusion.h"
#include "antioch/transport_mixture.h"
#include "antioch/physical_set.h"

// C++
#include <iostream>
#include <vector>


namespace Antioch
{
  // forward declaration
  template <typename CoeffType>
  class TransportSpecies;

  template<class NumericType, typename Interpolator>
  void build_molecular_binary_diffusion( PhysicalSet<MolecularBinaryDiffusion<NumericType, Interpolator> , TransportMixture<NumericType> >& D);

// ----------------------------------------- //

  template<class NumericType, class Interpolator>
  void build_molecular_binary_diffusion( PhysicalSet<MolecularBinaryDiffusion<NumericType, Interpolatir>, TransportMixture<NumericType> >& D)
  {
     Initializer<MolecularBinaryDiffusion<Numeric> > init;
     for(unsigned int s = 0; s < D.mixture().n_species(); s++)
     {
        init.si = D.mixture().transport_species()[s];
        for (unsigned int j = 0; j <= s; j++)
        {
           init.j = j;
           init.sj = D.mixture().transport_species()[j];
           D.add_model(D.mixture().species_inverse_name_map().at(s),init);
        }
     }
  }

}

#endif
