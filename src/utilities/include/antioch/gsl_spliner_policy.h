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

#ifndef ANTIOCH_GSL_SPLINER_POLICY_H
#define ANTIOCH_GSL_SPLINER_POLICY_H

// Antioch
#include "antioch/gsl_spliner_shim.h"

namespace Antioch
{
  namespace AntiochPrivate
  {

    template<bool B>
    struct GSLSplinerPolicy
    {
      template <typename Scalar>
      Scalar interpolation(const Scalar & x, const GSLSplinerShim& gsl_shim)
      {return gsl_shim.eval(x);}

      template <typename Scalar>
      Scalar dinterpolation(const Scalar & x, const GSLSplinerShim& gsl_shim)
      {return gsl_shim.eval_deriv(x);}
    };

    template <>
    struct GSLSplinerPolicy<true>
    {
      template <typename VectorScalar>
      VectorScalar interpolation(const VectorScalar & x, const GSLSplinerShim& gsl_shim)
      {
        VectorScalar out = zero_clone(x);
        for(unsigned int i =0; i < x.size(); ++i)
          {
            out[i] = gsl_shim.eval(x[i]);
          }
        return out;
      }

      template <typename VectorScalar>
      VectorScalar dinterpolation(const VectorScalar & x, const GSLSplinerShim& gsl_shim)
      {
        VectorScalar out = zero_clone(x);
        for(unsigned int i =0; i < x.size(); ++i)
          {
            out[i] = gsl_shim.eval_deriv(x[i]);
          }
        return out;
      }
    };

  } // end namespace AntiochPrivate

} // end namespace Antioch

#endif // ANTIOCH_GSL_SPLINER_POLICY_H
