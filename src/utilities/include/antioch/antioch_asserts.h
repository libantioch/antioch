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
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef ANTIOCH_ASSERTS_H
#define ANTIOCH_ASSERTS_H

// Antioch
#include "antioch/antioch_exceptions.h"
#include "antioch_config.h" // for ANTIOCH_HAVE_CXX11

// C++
#include <iostream>
#include <iomanip>

#define antioch_here()     do { std::cerr << __FILE__ << ", line " << __LINE__ << ", compiled " << __DATE__ << " at " << __TIME__ << std::endl; } while (0)

// The antioch_assert() macro acts like C's assert(), but throws a
// antioch_error() (including stack trace, etc) instead of just exiting

// When not debugging, we don't test any asserts
#ifdef NDEBUG
#define antioch_assert(asserted)  ((void) 0)
#define antioch_assert_msg(asserted, msg)  ((void) 0)
#define antioch_assert_equal_to(expr1,expr2)  ((void) 0)
#define antioch_assert_not_equal_to(expr1,expr2)  ((void) 0)
#define antioch_assert_less(expr1,expr2)  ((void) 0)
#define antioch_assert_greater(expr1,expr2)  ((void) 0)
#define antioch_assert_less_equal(expr1,expr2)  ((void) 0)
#define antioch_assert_greater_equal(expr1,expr2)  ((void) 0)
#else

#define antioch_assert(asserted)  do { if (!(asserted)) { std::cerr << "Assertion `" #asserted "' failed." << std::endl; antioch_error(); } } while(0)

// When using C++11, we can test asserts comparing two different types
// robustly
// #if __cplusplus > 199711L // http://gcc.gnu.org/bugzilla/show_bug.cgi?id=1773
#ifdef ANTIOCH_HAVE_CXX11
#define antioch_assert_equal_to(expr1,expr2)  do { typedef decltype(expr1) type1; typedef decltype(expr2) type2; if (!((expr1 == static_cast<type1>(expr2)) && static_cast<type2>(expr1) == expr2)) { std::cerr << "Assertion `" #expr1 " == " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; antioch_error(); } } while(0)
#define antioch_assert_not_equal_to(expr1,expr2)  do { typedef decltype(expr1) type1; typedef decltype(expr2) type2; if (!((expr1 != static_cast<type1>(expr2)) && (static_cast<type2>(expr1) != expr2))) { std::cerr << "Assertion `" #expr1 " != " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; antioch_error(); } } while(0)
#define antioch_assert_less(expr1,expr2)  do { typedef decltype(expr1) type1; typedef decltype(expr2) type2; if (!((static_cast<type2>(expr1) < expr2) && (expr1 < static_cast<type1>(expr2)))) { std::cerr << "Assertion `" #expr1 " < " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; antioch_error(); } } while(0)
#define antioch_assert_greater(expr1,expr2)  do { typedef decltype(expr1) type1; typedef decltype(expr2) type2; if (!((static_cast<type2>(expr1) > expr2) && (expr1 > static_cast<type1>(expr2)))) { std::cerr << "Assertion `" #expr1 " > " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; antioch_error(); } } while(0)
#define antioch_assert_less_equal(expr1,expr2)  do { typedef decltype(expr1) type1; typedef decltype(expr2) type2; if (!((static_cast<type2>(expr1) <= expr2) && (expr1 <= static_cast<type1>(expr2)))) { std::cerr << "Assertion `" #expr1 " <= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; antioch_error(); } } while(0)
#define antioch_assert_greater_equal(expr1,expr2)  do { typedef decltype(expr1) type1; typedef decltype(expr2) type2; if (!((static_cast<type2>(expr1) >= expr2) && (expr1 >= static_cast<type1>(expr2)))) { std::cerr << "Assertion `" #expr1 " >= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; antioch_error(); } } while(0)

// When using C++98, we let the compiler pick the type conversion and
// hope for the best.
#else
#define antioch_assert_equal_to(expr1,expr2)  do { if (!(expr1 == expr2)) { std::cerr << "Assertion `" #expr1 " == " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; antioch_error(); } } while(0)
#define antioch_assert_not_equal_to(expr1,expr2)  do { if (!(expr1 != expr2)) { std::cerr << "Assertion `" #expr1 " != " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; antioch_error(); } } while(0)
#define antioch_assert_less(expr1,expr2)  do { if (!(expr1 < expr2)) { std::cerr << "Assertion `" #expr1 " < " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; antioch_error(); } } while(0)
#define antioch_assert_greater(expr1,expr2)  do { if (!(expr1 > expr2)) { std::cerr << "Assertion `" #expr1 " > " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; antioch_error(); } } while(0)
#define antioch_assert_less_equal(expr1,expr2)  do { if (!(expr1 <= expr2)) { std::cerr << "Assertion `" #expr1 " <= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; antioch_error(); } } while(0)
#define antioch_assert_greater_equal(expr1,expr2)  do { if (!(expr1 >= expr2)) { std::cerr << "Assertion `" #expr1 " >= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; antioch_error(); } } while(0)

#endif // C++11

#endif // NDEBUG


#define antioch_error()    do { antioch_here(); ANTIOCH_THROW(Antioch::LogicError()); } while(0)
#define antioch_not_implemented()    do { antioch_here(); ANTIOCH_THROW(Antioch::NotImplemented()); } while(0)
#define antioch_file_error(filename)    do { antioch_here(); ANTIOCH_THROW(Antioch::FileError(filename)); } while(0)
#define antioch_unit_error(description)    do { antioch_here(); ANTIOCH_THROW(Antioch::UnitError(description)); } while(0)
#define antioch_parameter_required(parameter,defaultpar) \
         do { std::cerr <<  "\n*** Warning, The parameter " << parameter << " is not provided\n" \
                        <<  "default parameter is " << defaultpar << std::endl \
                        << __FILE__ << ", line " << __LINE__ << std::endl;} while(0)
#define antioch_unit_required(parameter,defaultunit) \
         do { std::cerr <<  "\n*** Warning, The parameter " << parameter << " is not given a unit\n" \
                        <<  "default unit given is " << defaultunit << std::endl \
                        << __FILE__ << ", line " << __LINE__ << std::endl;} while(0)
// The antioch_deprecated macro warns that you are using obsoleted code,
// savagely copy/pasted from libmesh
#define antioch_deprecated() \
          do  {std::cerr << "\n*** Warning, This code is deprecated, and likely to be removed in future library versions!\n" \
                         << __FILE__ << ", line " << __LINE__ << ", compiled " << __DATE__ << " at " << __TIME__ << " ***" << std::endl;} while(0)

#endif // ANTIOCH_ASSERTS_H
