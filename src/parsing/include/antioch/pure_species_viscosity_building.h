//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_PURE_SPECIES_BUILDING_H
#define ANTIOCH_PURE_SPECIES_BUILDING_H

// Antioch
#include "antioch/pure_species_viscosity.h"
#include "antioch/transport_mixture.h"
#include "antioch/physical_set.h"

// C++
#include <iostream>
#include <vector>


namespace Antioch
{
  template<class NumericType>
  void build_pure_species_viscosity( PhysicalSet<PureSpeciesViscosity<NumericType>, TransportMixture<NumericType> >& mu)

// ----------------------------------------- //

  template<class NumericType>
  void build_pure_species_viscosity( PhysicalSet<PureSpeciesViscosity<NumericType>, TransportMixture<NumericType> >& mu)
  {
      for(unsigned int s = 0; s < mu.mixture().n_species(); s++)
      {
          std::vector<CoeffType> coeffs(4,0);
          coeffs[0] = mu.mixture().transport_species()[s].LJ_depth();
          coeffs[1] = mu.mixture().transport_species()[s].LJ_diameter();
          coeffs[2] = mu.mixture().transport_species()[s].dipole_moment();
          coeffs[3] = mu.mixture().transport_species()[s].mass() / Constants::Avogadro<NumericType>();
          mu.add_model(mu.mixture().chemical_mixture().species_inverse_name_map().at(s),coeffs);
      }
  }

}

#endif
