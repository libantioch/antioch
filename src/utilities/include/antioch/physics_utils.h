//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-
//--------------------------------------------------------------------------

#ifndef ANTIOCH_PHYSICS_UTILS_H
#define ANTIOCH_PHYSICS_UTILS_H

// Antioch
#include "antioch/Stockmayer_potential.h"
#include "antioch/transport_mixture.h"

//C++

namespace Antioch
{
  //! check if the highest reduced temperature is within boundaries
  template <typename ThermoEvaluator,typename CoeffType>
  bool verify_highest_reduced_temperature(const TransportMixture<ThermoEvaluator, CoeffType> & mix, const StateType & T, StateType & Tmax);

  template <typename PhysicsModel, typename StateType>
  void extrapolate_Stockmayer(PhysicsModel & to_extrapolate, const StateType & T);

//---------------------------------

  template <typename ThermoEvaluator,typename CoeffType>
  inline
  bool verify_highest_reduced_temperature(const TransportMixture<ThermoEvaluator, CoeffType> & mix, const StateType & T, StateType & Tmax)
  {
     Tmax = Antioch::zero_clone(T);
     for(unsigned int s = 0; s < mix.n_species(); s++)
     {
        const TransportSpecies<CoeffType> & spec = *mix.transport_species()[s]; // convenient method

        if( Tmax < T / spec.LJ_depth() )Tmax = T/spec.LJ_depth();
     }

     return (Tmax > Stockmayer<CoeffType>().max_reduced_temperature());
  }

  template <typename PhysicsModel, typename StateType>
  inline
  void extrapolate_Stockmayer(PhysicsModel & to_extrapolate, const StateType & T)
  {
     to_extrapolate.built_extrapolation(T);
  }

}//end namespace Antioch


#endif //ANTIOCH_PHYSICS_UTILS_H
