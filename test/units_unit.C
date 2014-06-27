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

/* test on Units management
 *   - is_homogeneous
 *   - factor_to_some_unit
 *   - developSymbol
 *   - harmonizedSymbol
 *   - operations:
 *      - addition
 *      - substraction
 *      - multiplication
 *      - multiple combinations with parenthesises
 */

// C++
#include <limits>
#include <vector>

//Antioch
#include "antioch/units.h"
#include "antioch/unit_defs.h"

template <typename T>
int test_homogeneity()
{
  Antioch::Units<T> test("W");

  if(!test.is_homogeneous("J/s"))return 1;
  if(!test.is_homogeneous("N.m.s-1"))return 2;
  if(!test.is_homogeneous("kg.m2/s3"))return 3;
  if(!test.is_homogeneous("m2/kg-1.s-3"))return 4;
  if(!test.is_homogeneous("m.m.kg/s/s/s"))return 5;
  if(!test.is_homogeneous("erg/s"))return 6;

  std::cout << "Homogeneity done" << std::endl;

  return 0;

}

template <typename T>
int test_factor(int nAddition)
{
  using std::abs;


  const T tol = std::numeric_limits<T>::epsilon() * 2.;
  const T unity = 1.L;

// testing SI derived
  Antioch::Units<T> test("W");
  if(abs((test.get_SI_factor()-unity)/unity) > tol)
  {
    std::cerr << std::scientific << std::setprecision(16)
              << "mismatch in unit SI description, unity is not found" << std::endl
              << "factor = " << test.get_SI_factor() << std::endl;
    return 1;
  }
  test.set_unit("J");
  if(abs((test.get_SI_factor()-unity)/unity) > tol)
  {
    std::cerr << std::scientific << std::setprecision(16)
              << "mismatch in unit SI description, unity is not found" << std::endl
              << "factor = " << test.get_SI_factor() << std::endl;
    return 1;
  }
  test.set_unit("Hz");
  if(abs((test.get_SI_factor()-unity)/unity) > tol)
  {
    std::cerr << std::scientific << std::setprecision(16)
              << "mismatch in unit SI description, unity is not found" << std::endl
              << "factor = " << test.get_SI_factor() << std::endl;
    return 1;
  }
  test.set_unit("N");
  if(abs((test.get_SI_factor()-unity)/unity) > tol)
  {
    std::cerr << std::scientific << std::setprecision(16)
              << "mismatch in unit SI description, unity is not found" << std::endl
              << "factor = " << test.get_SI_factor() << std::endl;
    return 1;
  }
  test.set_unit("Pa");
  if(abs((test.get_SI_factor()-unity)/unity) > tol)
  {
    std::cerr << std::scientific << std::setprecision(16)
              << "mismatch in unit SI description, unity is not found" << std::endl
              << "factor = " << test.get_SI_factor() << std::endl;
    return 1;
  }
  test.set_unit("C");
  if(abs((test.get_SI_factor()-unity)/unity) > tol)
  {
    std::cerr << std::scientific << std::setprecision(16)
              << "mismatch in unit SI description, unity is not found" << std::endl
              << "factor = " << test.get_SI_factor() << std::endl;
    return 1;
  }
  test.set_unit("Bq");
  if(abs((test.get_SI_factor()-unity)/unity) > tol)
  {
    std::cerr << std::scientific << std::setprecision(16)
              << "mismatch in unit SI description, unity is not found" << std::endl
              << "factor = " << test.get_SI_factor() << std::endl;
    return 1;
  }

//testing composed unit
  srand(time(0)); //seed
  std::string uu("m");
  std::vector<std::string> add;
  for(int i = 0; i < nAddition; i++)
  {
    int indice = rand() % Antioch::UnitBaseStorage::known_units().n_known_units();
    add.push_back(Antioch::UnitBaseStorage::known_units().stored(indice).symbol());
  }

  int ndot(0);
  for(int i = 0; i < nAddition; i++)
  {
    int indice(0);
    if(ndot > 0)indice = rand() % ndot;
    size_t where(0);
    for(int j = 0; j < indice; j++)
    {
      where = uu.find(".",where);
    }
    if(where == 0)where = uu.size();
    std::string ad = "." + add[i];
    uu.insert(where,ad);

    if(ndot > 0)indice = rand() % ndot;
    where = 0;
    for(int j = 0; j < indice; j++)
    {
      where = uu.find(".",where);
    }
    if(where == 0)where = uu.size();
    ad = "/" + add[i];
    uu.insert(where,ad);
    ndot++;
  }

  test.set_unit(uu);
  if(abs((test.get_SI_factor()-unity)/unity) > tol)
  {
    std::cerr << std::scientific << std::setprecision(16)
              << "mismatch in unit combination, unity is not found" << std::endl
              << "factor = " << test.get_SI_factor() << std::endl;
    return 1;
  }
  
  std::cout << "Combination done at ";
  (nAddition < 100)?std::cout << "low":std::cout << "high"; 
  std::cout << " level" << std::endl;

  T RJ = 8.3144621L;
  T Rcal = 1.9858775L;
  T Runc = 0.0000034L;
  test.set_unit("J/mol/K");
  Antioch::Units<T> sec("cal/mol/K");
  T Rcalc = RJ * test.factor_to_some_unit("cal/mol/K");
  if(abs((Rcalc - Rcal)/Rcal) > tol && abs((Rcalc - Rcal)/Rcal) > abs(Runc/Rcal))
  {
    std::cerr << std::scientific << std::setprecision(16)
              << "mismatch in R calculation" << std::endl
              << "Rtheo = " << Rcal  << " cal/mol/K" << std::endl
              << "Rcal  = " << Rcalc << " cal/mol/K" << std::endl
              << "relative error  = " << abs((Rcalc - Rcal)/Rcal) << std::endl
              << "uncertainty in calorie = " << Runc << std::endl
              << "tol  = " << tol << std::endl;
    return 1;
  }

  std::cout << "R calculation passed" << std::endl;

  // now combining with parenthesises E = m * c * c, [E] = [m] + [c] + [c]
  Antioch::Units<T> m("g"), c("m/s");
  Antioch::Units<T> E = m + c + c;

  if(!E.is_homogeneous("J"))
  {
     std::cerr << "E = m * c *c failed, output unit is " << E.get_symbol() << std::endl;
     return 1;
  }

  std::cout << "E = m * c * c calculation passed" << std::endl;

  return 0;
}

template <typename T>
int tester()
{
  return (
          test_homogeneity<T>() ||
          test_factor<T>(10)    ||
          test_factor<T>(200)
          );
}


int main()
{
  return (
          tester<float>() ||
          tester<double>()    ||
          tester<long double>()
          );
}
