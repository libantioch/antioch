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

// C++
#include <iostream>
#include <iomanip>

#define antioch_here()     do { std::cerr << __FILE__ << ", line " << __LINE__ << ", compiled " << __DATE__ << " at " << __TIME__ << std::endl; } while (0)

// The antioch_assert() macro acts like C's assert(), but throws a
// antioch_error() (including stack trace, etc) instead of just exiting
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
#define antioch_assert_equal_to(expr1,expr2)  do { if (!(expr1 == expr2)) { std::cerr << "Assertion `" #expr1 " == " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; antioch_error(); } } while(0)
#define antioch_assert_not_equal_to(expr1,expr2)  do { if (!(expr1 != expr2)) { std::cerr << "Assertion `" #expr1 " != " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; antioch_error(); } } while(0)
#define antioch_assert_less(expr1,expr2)  do { if (!(expr1 < expr2)) { std::cerr << "Assertion `" #expr1 " < " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; antioch_error(); } } while(0)
#define antioch_assert_greater(expr1,expr2)  do { if (!(expr1 > expr2)) { std::cerr << "Assertion `" #expr1 " > " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; antioch_error(); } } while(0)
#define antioch_assert_less_equal(expr1,expr2)  do { if (!(expr1 <= expr2)) { std::cerr << "Assertion `" #expr1 " <= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; antioch_error(); } } while(0)
#define antioch_assert_greater_equal(expr1,expr2)  do { if (!(expr1 >= expr2)) { std::cerr << "Assertion `" #expr1 " >= " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; antioch_error(); } } while(0)
#endif


#define antioch_error()    do { antioch_here(); ANTIOCH_THROW(Antioch::LogicError()); } while(0)
#define antioch_not_implemented()    do { antioch_here(); ANTIOCH_THROW(Antioch::NotImplemented()); } while(0)
#define antioch_file_error(filename)    do { antioch_here(); ANTIOCH_THROW(Antioch::FileError(filename)); } while(0)

#endif // ANTIOCH_ASSERTS_H
