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
//
// $Id: valarray_utils.h 37170 2013-02-19 21:40:39Z roystgnr $
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef ANTIOCH_VEXCL_UTILS_DECL_H
#define ANTIOCH_VEXCL_UTILS_DECL_H

#ifdef ANTIOCH_METAPROGRAMMING_H
#  error valarray_utils_decl.h must be included before metaprogramming.h
#endif

// Antioch
#include "antioch/metaprogramming_decl.h"

#include <type_traits> // std::enable_if

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
  template <typename T> class vector_expression;
  template <typename T> class is_vector_expression;
}
namespace boost {
  namespace proto {
    namespace result_of {
      template <typename T1, typename T2, typename T3, typename T4, typename T5>
      class make_expr;
    }
    namespace tag {
      class if_else_;
    }
    class deduce_domain;
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
ANTIOCH_RETURNEXPR(vex::vector<T>, vex::max(a,b));

template <typename T>
inline
ANTIOCH_AUTO(vex::vector<T>)
min(const vex::vector<T>& a,
    const vex::vector<T>& b)
ANTIOCH_RETURNEXPR(vex::vector<T>, vex::min(a,b));
}
#endif


namespace Antioch
{

// Class to allow tag dispatching to VexCL specializations
struct vexcl_library_tag : public numeric_library_tag {};

template <typename T>
struct tag_type<T,
          typename std::enable_if<vex::is_vector_expression<T>::value>::type
        >
{
    typedef const vexcl_library_tag type;
};

template <typename T, typename NewScalar>
struct rebind<vex::vector<T>, NewScalar>
{
  typedef vex::vector<NewScalar> type;
};

template <typename T>
inline
T
max (const vex::vector<T>& in);

template <typename T>
inline
T
min (const vex::vector<T>& in);

template <typename T>
struct return_auto<vex::vector<T> >;

template <typename T>
struct has_size<vex::vector<T> >;

template <typename T>
struct size_type<vex::vector<T> >;

template <typename T>
struct value_type<vex::vector<T> >;

template <typename T>
struct raw_value_type<vex::vector<T> >;

template <typename T>
inline
vex::vector<T>
zero_clone(const vex::vector<T>& example);

template <typename T1, typename T2>
inline
void
zero_clone(vex::vector<T1>& output, const vex::vector<T2>& example);

template <typename T, typename Scalar>
inline
vex::vector<T>
constant_clone(const vex::vector<T>& example, const Scalar& value);

template <typename T>
inline
void
init_clone(vex::vector<T>& output, const vex::vector<T>& example);

template <typename BoolInput, typename IfValue, typename ElseValue>
typename boost::proto::result_of::make_expr<
  boost::proto::tag::if_else_,
  boost::proto::deduce_domain,
  const vex::vector_expression<BoolInput>&,
  const IfValue&,
  const ElseValue&
>::type const
if_else(const vex::vector_expression<BoolInput> &condition,
	const IfValue   &if_true,
	const ElseValue &if_false);

template <typename T>
inline
bool disjunction_root(const T & vec_input, vexcl_library_tag);

template <typename T>
inline
bool conjunction_root(const T & vec_input, vexcl_library_tag);


template <typename VectorT, typename IntT>
inline
typename enable_if_c<
  vex::is_vector_expression<typename value_type<VectorT>::type>::value &&
  vex::is_vector_expression<IntT>::value,
  typename value_type<VectorT>::type
>::type
eval_index(const VectorT& vec, const IntT& index);

} // end namespace Antioch


#endif //ANTIOCH_VEXCL_UTILS_DECL_H
