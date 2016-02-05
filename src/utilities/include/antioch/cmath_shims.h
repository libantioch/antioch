
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014-2016 Paul T. Bauman, Benjamin S. Kirk,
//                         Sylvain Plessis, Roy H. Stonger
//
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

#ifndef ANTIOCH_CMATH_H
#define ANTIOCH_CMATH_H

// Antioch headers
#include "antioch/metaprogramming_decl.h"

// C++ headers
#include <algorithm> // max, min
#include <cmath>     // everything else

// Create shim methods Antioch::ant_exp() for exp(), etcetera.

namespace Antioch
{
  // Bring std:: methods solely into the shim methods.
  // There has to be a better way to trick an auto function into doing
  // Koenig lookup on its return type...

#ifdef ANTIOCH_HAVE_CXX11
#define ANTIOCH_UNARY_SHIM(funcname) \
namespace Ugly_CXX11_Workaround \
{ \
  using std::funcname; \
  template <typename T> \
  inline \
  auto \
  ant_##funcname (const T& in) \
  -> decltype (funcname(in)) \
  { return funcname(in); } \
} \
using Ugly_CXX11_Workaround::ant_##funcname;

#define ANTIOCH_BINARY_SHIM(funcname) \
namespace Ugly_CXX11_Workaround \
{ \
  using std::funcname; \
  template <typename T1, typename T2> \
  inline \
  auto \
  ant_##funcname (const T1& in1, const T2& in2) \
  -> decltype (funcname(in1, in2)) \
  { return funcname(in1, in2); } \
} \
using Ugly_CXX11_Workaround::ant_##funcname;

#else // ANTIOCH_HAVE_CXX11

#define ANTIOCH_UNARY_SHIM(funcname) \
  template <typename T> \
  inline \
  T \
  ant_##funcname (const T& in) \
  { \
    using std::funcname; \
    return funcname(in); \
  }

#define ANTIOCH_BINARY_SHIM(funcname) \
  template <typename T1, typename T2> \
  inline \
  T1 \
  ant_##funcname (const T1& in1, const T2& in2) \
  { \
    using std::funcname; \
    return funcname(in1, in2); \
  }

#endif // ANTIOCH_HAVE_CXX11

ANTIOCH_UNARY_SHIM(exp)
ANTIOCH_UNARY_SHIM(log)
ANTIOCH_UNARY_SHIM(log10)
ANTIOCH_UNARY_SHIM(sin)
ANTIOCH_UNARY_SHIM(cos)
ANTIOCH_UNARY_SHIM(tan)
ANTIOCH_UNARY_SHIM(asin)
ANTIOCH_UNARY_SHIM(acos)
ANTIOCH_UNARY_SHIM(atan)
ANTIOCH_UNARY_SHIM(sinh)
ANTIOCH_UNARY_SHIM(cosh)
ANTIOCH_UNARY_SHIM(tanh)
ANTIOCH_UNARY_SHIM(sqrt)
ANTIOCH_UNARY_SHIM(abs)
ANTIOCH_UNARY_SHIM(fabs)
ANTIOCH_UNARY_SHIM(ceil)
ANTIOCH_UNARY_SHIM(floor)

ANTIOCH_BINARY_SHIM(pow)
ANTIOCH_BINARY_SHIM(atan2)
ANTIOCH_BINARY_SHIM(max)
ANTIOCH_BINARY_SHIM(min)
ANTIOCH_BINARY_SHIM(fmod)

} // end namespace Antioch

#endif //ANTIOCH_CMATH_H
