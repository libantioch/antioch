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
 *   - isHomogeneous
 *   - FactorToSomeUnit
 *   - developSymbol
 *   - harmonizedSymbol
 *   - operations:
 *      - addition
 *      - substraction
 *      - multiplication
 */

#include "antioch/Units.hpp"
#include "antioch/unit_defs.hpp"
#include <vector>

int test_homogeneity()
{
  Antioch::Units test("W");

  if(!test.isHomogeneous("J/s"))return 1;
  if(!test.isHomogeneous("N.m.s-1"))return 2;
  if(!test.isHomogeneous("kg.m2/s3"))return 3;
  if(!test.isHomogeneous("m2/kg-1.s-3"))return 4;
  if(!test.isHomogeneous("m.m.kg/s/s/s"))return 5;
  if(!test.isHomogeneous("erg/s"))return 6;

  std::cout << "Homogeneity done" << std::endl;

  return 0;

}

int test_factor(int nAddition)
{
// testing SI derived
  Antioch::Units test("W");
  if(test.getSIFactor() != 1.)return 1;
  test.setUnit("J");
  if(test.getSIFactor() != 1.)return 2;
  test.setUnit("Hz");
  if(test.getSIFactor() != 1.)return 3;
  test.setUnit("N");
  if(test.getSIFactor() != 1.)return 4;
  test.setUnit("Pa");
  if(test.getSIFactor() != 1.)return 5;
  test.setUnit("C");
  if(test.getSIFactor() != 1.)return 6;
  test.setUnit("Bq");
  if(test.getSIFactor() != 1.)return 7;

//testing composed unit
  srand(time(0)); //seed
  std::string uu("m");
  std::vector<std::string> add;
  for(int i = 0; i < nAddition; i++)
  {
    int indice = rand() % Antioch::nKnownUnits;
    add.push_back(Antioch::knownUnits[indice].getSymbol());
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

  test.setUnit(uu);
  if((test.getSIFactor() - 1.) > 1e-15)return 8;
  
  std::cout << "Combination done at ";
  (nAddition < 100)?std::cout << "low":std::cout << "high"; 
  std::cout << " level" << std::endl;

  return 0;
}

int main()
{
  return (
          test_homogeneity() ||
          test_factor(10)    ||
          test_factor(200)
          );
}
