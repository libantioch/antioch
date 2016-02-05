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

#include "antioch_config.h"

#ifdef ANTIOCH_HAVE_CPPUNIT
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
#endif // ANTIOCH_HAVE_CPPUNIT

int main()
{
#ifdef ANTIOCH_HAVE_VEXCL
#ifdef ANTIOCH_HAVE_CPPUNIT
  CppUnit::TextUi::TestRunner runner;
  CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
  runner.addTest( registry.makeTest() );

  // If the tests all succeed, report success
  if (runner.run())
    return 0;

  // If any test fails report failure
  return 1;
#else // ANTIOCH_HAVE_CPPUNIT
  // If we don't have CPPUnit, report we skipped
  // 77 return code tells Automake we skipped this.
  return 77;
#endif // ANTIOCH_HAVE_CPPUNIT

#else // ANTIOCH_HAVE_VEXCL
  // If we don't have VexCL, report we skipped
  // 77 return code tells Automake we skipped this.
  return 77;
#endif // ANTIOCH_HAVE_VEXCL
}
