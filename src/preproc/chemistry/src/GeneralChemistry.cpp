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
#include "GeneralChemistry.hpp"
#include "coreObjects.hpp"
#include <sstream>

std::vector<std::string> splitInAtoms(const std::string &molecule)
{
  const std::string method("std::vector<std::string> splitInAtoms(const std::string &)");

  std::vector<std::string> atoms;
  std::string atom,num;
  for(unsigned int i = 0; i < molecule.size(); i++)
  {
      char c = molecule[i];
      if(c == '0' || c == '1' || c == '2' ||
         c == '3' || c == '4' || c == '5' ||
         c == '6' || c == '7' || c == '8' ||
         c == '9') //if figure
      {
        if(atom.empty())continue;
        num.push_back(c);        
      }else //letter
      {
        if(!num.empty())
        {
          if(atom.empty())CoreError err("Error",method,"This error is not supposed to happen, ever, never, don't you think of it!!");
          std::stringstream nat(num);
          int natom;
          nat >> natom;
          for(int i = 0; i < natom; i++)
          {
             atoms.push_back(atom);
          }
          num.clear();
          atom.clear();
        }else if(c == 'A' || c == 'B' || c == 'C' || c == 'D' || c == 'E' || 
                 c == 'F' || c == 'G' || c == 'H' || c == 'I' || c == 'J' || 
                 c == 'K' || c == 'L' || c == 'M' || c == 'N' || c == 'O' || 
                 c == 'P' || c == 'Q' || c == 'R' || c == 'S' || c == 'T' || 
                 c == 'U' || c == 'V' || c == 'W' || c == 'X' || c == 'Y' || c == 'Z')
        {
          if(!atom.empty())
          {
             atoms.push_back(atom);
             atom.clear();
          }
        }else if(c == '(')//state
        {
          while(molecule[i] != ')')
          {
             i++;
             if(i >= molecule.size())CoreError err("Error",method,"An open parenthesis is not matched in your molecule " + molecule);
          }
          continue;
        }else if(c == '[')//isotopomere
        {
          while(molecule[i] != ']')
          {
             i++;
             if(i >= molecule.size())CoreError err("Error",method,"An open square braquet is not matched in your molecule " + molecule);
          }
          continue;
        }
        atom.push_back(c);
      }
  }
  if(atom.empty())CoreError err("Error",method,"This molecule has a strange name (" + molecule + ")");
  int natom(1);
  if(!num.empty())
  {
    std::stringstream nat(num);
    nat >> natom;
  }
  for(int i = 0; i < natom; i++)
  {
    atoms.push_back(atom);
  }
  num.clear();
  atom.clear();

  return atoms;
}

bool isStrictIsomere(const std::string &mol1, const std::string &mol2)
{
  std::vector<std::string> atoms1 = splitInAtoms(mol1);
  std::vector<std::string> atoms2 = splitInAtoms(mol2);

  if(atoms1.size() != atoms2.size())return false;

  for(unsigned int i = 0; i < atoms1.size(); i++)
  {
     if(atoms1[i] != atoms2[i])return false;
  }

  return true;
}


bool isIsomere(const std::string &mol1, const std::string &mol2)
{
  std::vector<std::string> atoms1 = splitInAtoms(mol1);
  std::vector<std::string> atoms2 = splitInAtoms(mol2);

  if(atoms1.size() != atoms2.size())return false;

  unsigned int nj(atoms2.size());
  for(unsigned int i = 0; i < atoms1.size(); i++)
  {
     for(unsigned int j = 0; j < atoms2.size(); j++)
     {
        if(atoms1[i] == atoms2[j])
        {
           atoms2.erase(atoms2.begin() + j);
        }
     }
     if(atoms2.size() == nj)return false;
     nj--;
  }
  return true;
}
