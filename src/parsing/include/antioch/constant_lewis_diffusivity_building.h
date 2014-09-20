//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_CONSTANT_LEWIS_DIFFUSIVITY_BUILDING_H
#define ANTIOCH_CONSTANT_LEWIS_DIFFUSIVITY_BUILDING_H

// Antioch
#include "antioch/constant_lewis_diffusivity.h"

// C++
#include <iostream>
#include <vector>


namespace Antioch
{

  //
  template <typename NumericType>
  class ChemicalMixture;


  template<class NumericType>
  void build_constant_lewis_diffusivity( PhysicalSet<ConstantLewisDiffusivity<NumericType> , ChemicalMixture<NumericType> >& D, const  NumericType & Le);

// ----------------------------------------- //

  template<class NumericType>
  void build_constant_lewis_diffusivity( PhysicalSet<ConstantLewisDiffusivity<NumericType>, ChemicalMixture<NumericType> >& D, const NumericType & Le)
  {
          D.set() = new ConstantLewisDiffusivity<NumericType>(Le);
  }

}

#endif