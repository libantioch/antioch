//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
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
//
// $Id: valarray_utils.h 37170 2013-02-19 21:40:39Z roystgnr $
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef ANTIOCH_VEXCL_UTILS_H
#define ANTIOCH_VEXCL_UTILS_H

#ifdef ANTIOCH_METAPROGRAMMING_H
#  ifndef ANTIOCH_VEXCL_UTILS_DECL_H
#    error valarray_utils_decl.h must be included before metaprogramming.h
#  endif
#endif

// Antioch
#include "antioch/metaprogramming.h"
#include "antioch/antioch_asserts.h"


#ifdef ANTIOCH_HAVE_VEXCL
// Though the following implementations are all valid without
// <vexcl/vexcl.hpp>, successfully using them with
// VexCL types requires VexCL to be included first.
// Configure-time VexCL support enforces this constraint but
// header-only VexCL may be mixed with header-only Antioch
// without configure-time flags.
#include "vexcl/vexcl.hpp"
#else
// Forward declaration instead
namespace vex {
template <typename T> class vector;

template<typename real, class RDC>
class Reductor;

class MAX;
class MIN;
}
namespace boost {
  namespace proto {
    template<typename A0, typename A1, typename A2> 
      typename result_of::make_expr<
        tag::if_else_, 
        deduce_domain,
        A0 const &, 
        A1 const &, 
        A2 const &
      >::type const 
      if_else(A0 const & a0, A1 const & a1, A2 const & a2);
  }
}
#endif


#ifdef VEXCL_OPERATIONS_HPP
namespace std {
template <typename T>
inline
ANTIOCH_AUTO(vex::vector<T>)
max(const vex::vector<T>& a,
    const vex::vector<T>& b)
ANTIOCH_AUTOFUNC(vex::vector<T>, vex::max(a,b))

template <typename T>
inline
ANTIOCH_AUTO(vex::vector<T>)
min(const vex::vector<T>& a,
    const vex::vector<T>& b)
ANTIOCH_AUTOFUNC(vex::vector<T>, vex::min(a,b))
}

#endif


namespace Antioch
{

template <typename T>
inline
T
max (const vex::vector<T>& in)
{
  vex::Reductor<double, vex::MAX> vex_max(in.queue_list());
  return vex_max(in);
}

template <typename T>
inline
T
min (const vex::vector<T>& in)
{
  vex::Reductor<double, vex::MIN> vex_min(in.queue_list());
  return vex_min(in);
}

template <typename T>
struct return_auto<vex::vector<T> >
{
  static const bool value = true;
};

template <typename T>
struct has_size<vex::vector<T> >
{
  static const bool value = true;
};

template <typename T>
struct size_type<vex::vector<T> >
{
  typedef std::size_t type;
};

template <typename T>
struct value_type<vex::vector<T> >
{
  typedef T type;
};

template <typename T>
struct raw_value_type<vex::vector<T> >
{
  typedef typename raw_value_type<T>::type type;
};

template <typename T>
inline
vex::vector<T>
zero_clone(const vex::vector<T>& example)
{
  vex::vector<T> returnval(example.queue_list(), example.size());
  returnval.clear();
  return returnval;
}

template <typename T1, typename T2>
inline
void
zero_clone(vex::vector<T1>& output, const vex::vector<T2>& example)
{
  output.resize(example.queue_list(), example.size());
  output = 0;
}

template <typename T, typename Scalar>
inline
vex::vector<T>
constant_clone(const vex::vector<T>& example, const Scalar& value)
{
  vex::vector<T> returnval(example.queue_list(), example.size());
  returnval = value;
  return returnval;
}

template <typename T>
inline
void
init_clone(vex::vector<T>& output, const vex::vector<T>& example)
{
  output.resize(example.queue_list(), example.size());
  output = example;
}

template <typename BoolInput,
	  typename IfValue, typename ElseValue>
inline
typename boost::proto::result_of::make_expr<
  boost::proto::tag::if_else_,
  boost::proto::deduce_domain,
  const vex::vector_expression<BoolInput>&,
  const IfValue&,
  const ElseValue&
>::type const
if_else(const vex::vector_expression<BoolInput> &condition,
	const IfValue   &if_true,
	const ElseValue &if_false)
{
  return boost::proto::if_else(condition, if_true, if_false);
}

} // end namespace Antioch


#endif //ANTIOCH_VEXCL_UTILS_H
