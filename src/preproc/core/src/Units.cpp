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
//Antioch
#include "antioch/Units.hpp"
#include "antioch/unit_defs.hpp"

//C++
#include <vector>
#include <cmath>

namespace Antioch{

Units::Units(std::string sym,std::string na,
             double conva,double convb,
             int mi,int kgi, int si, int Ai, int Ki, int moli, int cdi, int radi): //fully descriptive constructor
  symbol(sym),
  name(na),
  toSI(conva,convb),
  power(mi,kgi,si,Ai,Ki,moli,cdi,radi)
{
  developSymbol(symbol);
}

Units::Units(std::string sym,std::string na):  //constructors when automatic conversion
       symbol(sym),
       name(na)
{
  developSymbol(symbol);
  fillInPower(true);
}

Units::Units(std::string sym,Converter conv, std::string na)://constructor for given conversion
       symbol(sym),
       name(na),
       toSI(conv)
{
  developSymbol(symbol);
  fillInPower(false);
}

void Units::fillInPower(bool doConv)
{
  if(symbol.empty())return; // no unity
  const std::string method("Units::fillInPower()");
  std::string tmp(""),symboltmp(symbol);
  int signe(1),istart(0);

  while(symboltmp != contractedSymbol(symboltmp))symboltmp = contractedSymbol(symboltmp);
  if(symboltmp.empty())return; // no unity

  if(symboltmp[0] == '/')
  {
    signe = -1;
    istart = 1;
  }
  for(unsigned int i = istart ; i < symboltmp.size() ; i++)
  {
    if(symboltmp[i] == '.')
    {
      if(!parseSingleUnit(signe,tmp,doConv))
      {
        antiochError(method,"In symbol " + symboltmp + ", unit \"" + tmp + "\" does not ring a bell");
        break;
      }
      signe = 1;
      tmp.clear();
      }else if(symboltmp[i] == '/')
      {
       if(!parseSingleUnit(signe,tmp,doConv))
       {
         antiochError(method,"In symbol " + symboltmp + ", unit \"" + tmp + "\" does not ring a bell");
         break;
       }
       signe = -1;
       tmp.clear();
       }else
       {
         tmp += symboltmp[i];
       }
   }
   if(!parseSingleUnit(signe,tmp,doConv))
      antiochError(method,"In symbol " + symboltmp + ", unit \"" + tmp + "\" does not ring a bell");
}

void Units::parsePrefixeUnit(int &iUnit,int &iPre,std::string unit) const
{
   for(iPre = 0 ; iPre < nSIPrefixes ; iPre++)
   {
    if(unit.size() < Prefixes[iPre].symbol().size())continue;
    if(unit.substr(0,Prefixes[iPre].symbol().size()) != Prefixes[iPre].symbol())continue; //is it a unit ?

    iUnit = indexUnit(unit.substr(Prefixes[iPre].symbol().size(),std::string::npos)); //is it a unit ?

    if(unit.substr(0,Prefixes[iPre].symbol().size()) == Prefixes[iPre].symbol() && //if this is a prefixe AND
       iUnit != -1)
    {
      return; //this is a unit
    }
   }

   iPre = -1; //no prefixe
   iUnit = indexUnit(unit); //is it a unit ?

}

bool Units::parseSingleUnit(int signe,std::string unit,bool doConv)
{
  int iUnit = -1, nc = 0;
  bool iostat(false);
  int ipower = parsePower(unit,nc); // find power
  int iPre = -1;

  unit = unit.substr(0,unit.length()-nc); //unit without power

  parsePrefixeUnit(iUnit,iPre,unit);//find prefixe and unit

  if(iUnit == -1)return iostat; //not found

  double pre = 1.;
  if(iPre != -1)pre = Prefixes[iPre].value();
  if(knownUnits[iUnit].getSymbol() == "kg")pre *= 1e-3;
  InSI powerTmp = knownUnits[iUnit].getPower();
  power += (powerTmp * signe * ipower);
  if(doConv)
  {
    Converter convTmp = knownUnits[iUnit].getSIcoef() * pre;
    toSI += raise(convTmp,signe * ipower);
  }
  iostat = true;

  return iostat;
}

int Units::indexUnit(std::string unit) const
{
  int iUnit;

  if(unit.empty())return -1;

  if(unit == "g")//special case to adapt SI
  {
    unit = "kg";
  }
  for (iUnit = 0; iUnit < nKnownUnits ; iUnit++)
  {
    if(unit == knownUnits[iUnit].getSymbol())break;
  }

  if(iUnit >= nKnownUnits)iUnit=-1;


  return iUnit;
}

int Units::parsePower(std::string unit,int &nc) const
{
  const std::string method("Units::parsePower()");
  int ip = 1, loc = unit.length();
  char c = unit[loc-1];
  nc = 0;
  while(isNumber(c)){
    nc++;
    loc--;
    c = unit[loc-1];
  }

  if(unit[loc-1] == '-'){
        loc--;
        nc++;
  }

  std::string power = unit.substr(loc,std::string::npos);
  if(power.size() > 0)
  {
     std::stringstream p;
     p << power;
     p >> ip;
  }

  if(ip == 0)antiochError(method,"Invalid power found: " + unit);

  return ip;
}

const bool Units::isHomogeneous(std::string target) const
{
  if(target.empty())
  {
    return (!isUnited());
  }else
  {
    Units rhs(target);
    return isHomogeneous(rhs);
  }
  
}

int Units::getSIPower(const std::string &SIask) const
{
  const std::string method("int Units::getSIPower(const std::string &) const");
  if(SIask == "m")  return power.get_m();
  if(SIask == "kg") return power.get_kg();
  if(SIask == "s")  return power.get_s();
  if(SIask == "A")  return power.get_A();
  if(SIask == "K")  return power.get_K();
  if(SIask == "mol")return power.get_mol();
  if(SIask == "cd") return power.get_cd();
  if(SIask == "rad")return power.get_rad();

  antiochError(method,SIask + " is not a SI symbol");

  return 0;
}

std::string const Units::getSISymb() const
{
 std::string SISymb;
 SISymb.clear();

 if(power.get_m() != 0)SISymb += addSI(power.get_m(),"m");

 if(power.get_kg() != 0)
 {
   if(!SISymb.empty())SISymb += ".";
   SISymb += addSI(power.get_kg(),"kg");
 }
 if(power.get_s() != 0)
 {
   if(!SISymb.empty())SISymb += ".";
   SISymb += addSI(power.get_s(),"s");
 }
 if(power.get_A() != 0)
 {
   if(!SISymb.empty())SISymb += ".";
   SISymb += addSI(power.get_A(),"A");
 }
 if(power.get_K() != 0)
 {
   if(!SISymb.empty())SISymb += ".";
   SISymb += addSI(power.get_K(),"K");
 }
 if(power.get_mol() != 0)
 {
   if(!SISymb.empty())SISymb += ".";
   SISymb += addSI(power.get_mol(),"mol");
 }
 if(power.get_cd() != 0)
 {
   if(!SISymb.empty())SISymb += ".";
   SISymb += addSI(power.get_cd(),"cd");
 }
 if(power.get_rad() != 0)
 {
   if(!SISymb.empty())SISymb += ".";
   SISymb += addSI(power.get_rad(),"rad");
 }

 return SISymb;

}

const std::string Units::addSI(int pow,std::string symb) const
{
 std::string out("");

 if(pow == 0)return out;

 out = symb;
 if(pow != 1)
 {
   std::stringstream po;
   po << pow;
   out += po.str();
 }
 
 return out;

}

std::string const Units::manipulateSymbol(std::string input, bool contract) const
{
  const std::string method("Units::manipulateSymbol(const std::string &, bool) const");
  std::string harmSymb("");
  if(input.empty())input = symbol;

  std::vector<std::string> unitvec;
  std::vector<int> powervec;
  std::string curUnit(""),interUnit(".");//,strPower;

  for(unsigned int i = 0; i < input.size(); i++)
  {
    if(input[i] != '.' && input[i] != '/')curUnit += input[i];
    if(input[i] == '.' || input[i] == '/' || i == input.size() - 1 )
    {
      if(curUnit.empty())continue;
//first, parsing the unit
      int nc(0);
      int iUnit(-1),iPre(-1);
      int curPower = parsePower(curUnit,nc); //power
//      strPower = curUnit.substr(curUnit.length() - nc,std::string::npos); //resulting power
      curUnit = curUnit.substr(0,curUnit.length() - nc); //resulting unit
      parsePrefixeUnit(iUnit,iPre,curUnit);
 // unit is Prefixes[iPre].symbol() + knownUnits[iUnit].getSymbol() + strPower
      if(iUnit == -1)
      {
        antiochError(method,"The unit \"" + curUnit + "\" is unknown. No harmonized symbol will be produced.");
        harmSymb.clear();
        return harmSymb;
      }
//updating variables:
//  - curUnit is the unit (std::string)
//  - curPower is the power (int)
      if(interUnit == "/")curPower *= -1;

// checking what we've got
      unsigned int j;
      for(j = 0; j < unitvec.size(); j++)
      {
         Units tmp(unitvec[j]);
         bool same(tmp.getSymbol() == curUnit);//only contraction: strong condition
         if(!contract && !same)same = tmp.isHomogeneous(curUnit); //harmonizing: weak condition
         if(same)//if in there, update the power, 
         {
            powervec[j] += curPower;
            break;
         }
      }
      if(j >= unitvec.size())//if not, add
      {
        unitvec.push_back(curUnit);
        powervec.push_back(curPower);
      }
      interUnit = input[i];
      curUnit.clear();
    }
  }

//litre fix: only if harmonizing
//  if there is l AND m in the vector, all is converted to m. 
  if(!contract)
  {
    bool ism(false),isl(false);
    unsigned int herem(0),herel(0);
    for(unsigned int i = 0; i < unitvec.size(); i++)
    {
      if(unitvec[i] == "m")
      {
        ism = true;
        herem = i;
      }
      if(unitvec[i] == "l")
      {
        isl = true;
        herel = i;
      }
    }
    if(ism && isl)
    {
      powervec[herem] += 3 * powervec[herel];
      unitvec.erase(unitvec.begin() + herel);
      powervec.erase(powervec.begin() + herel);
    }
  }
//now writing harmsymb

  std::ostringstream outsym;
  int k(0);
  for(unsigned int i = 0; i < unitvec.size(); i++)
  {
    if(powervec[i] == 0)continue; //ignore deleted unit
    if(k != 0) // need a symbol, '.' or '/'
    {
      if(powervec[i] > 0) // '.'
      {
        outsym << ".";
      }else if(powervec[i] < 0) // '/' and reverse the power
      {
        outsym << "/";
        powervec[i] *= -1;
      }
    }
    outsym << unitvec[i];
    if(powervec[i] != 1)outsym << powervec[i];
    k++;
  }

  harmSymb = outsym.str();

  return harmSymb;
}

const int Units::nDimensionOfUnits() const
{
  int ndim(0);

  if(power.get_m()   != 0)ndim++;
  if(power.get_kg()  != 0)ndim++;
  if(power.get_s()   != 0)ndim++;
  if(power.get_A()   != 0)ndim++;
  if(power.get_K()   != 0)ndim++;
  if(power.get_mol() != 0)ndim++;
  if(power.get_cd()  != 0)ndim++;
  if(power.get_rad() != 0)ndim++;

  return ndim;
}

const std::string Units::getSIConvenientSymb() const
{
  if(nDimensionOfUnits() == 0)
  {
    std::string outStr("");
    return outStr;
  }else
  {
    for(int iUnit = 1; iUnit < nKnownUnits; iUnit++)
    {
      if(isHomogeneous(knownUnits[iUnit]) && 
         knownUnits[iUnit].getSIFactor() == 1. && 
         knownUnits[iUnit].getSITranslator() == 0.)return knownUnits[iUnit].getSymbol();
    }
    return getSISymb();
  }
}

void Units::setUnit(const std::string &sym, std::string na)
{
  symbol = sym;
  developSymbol(symbol);
  name = na;
  toSI.clear();
  power.clear();
  fillInPower(true);
}

void Units::developSymbol(std::string &subsymb)
{
  if(subsymb == "no unit" || subsymb == "No unit" || subsymb == "NO UNIT")
  {
     subsymb.clear();
     return;
  }
  if(subsymb.find("(") == std::string::npos)return;

//model is a serie ().()/()..., we find it
// we count the opening parenthesises (no) and
// the closing parenthesises (nc). We need to
// keep track of the position of opening (po)
// and closing (pc) parenthesises. We have a pair
// when no == nc.
  unsigned int no(0),nc(0);
  unsigned int po(0),pc(subsymb.size() - 1);
  for(unsigned int cc = 0; cc < subsymb.size(); cc++)
  {
    if(subsymb[cc] == '(')
    {
      no++;
      if(po == 0 && cc != 0)po = cc;
    }
    if(subsymb[cc] == ')')
    {
      nc++;
      pc = cc;
    }
    if(no == 0)continue;

    if(no == nc)//found a pair
    {
     //develop it
      if(pc == po + 1)
      {
        unsigned int off(0);
        if(po == 0)off = 1;//if ( is the first character or not
        subsymb.erase(po - 1 + off,3 - off);//if yes, suppress '()', if not suppress '.()' or '/()'
      }else
      {
        std::string insideStr = subsymb.substr(po + 1,pc - po - 1);
        developSymbol(insideStr);
        if(po != 0)if(subsymb[po - 1] == '/')reversePowerSymbol(insideStr);
        subsymb.replace(po,pc - po + 1,insideStr);
      }
     //reset the system
      no = 0;
      nc = 0;
      po = 0;
      pc = subsymb.size() - 1;
      cc -= 2;
    }
  }

//if first character is a power determinant ('/' or '.')
  if(subsymb[0] == '/') //we change only the first atomic unit
  {
    subsymb.erase(0,1);
    po = subsymb.find(".");
    if(po > subsymb.find("/"))po = subsymb.find("/");
    std::string curUnit = subsymb.substr(0,po);
    int nc(0),pow = - parsePower(curUnit,nc); //power
    if(pow != 1)
    {
      std::ostringstream np;
      np << pow;
      subsymb.replace(po - nc,nc,np.str());
    }else
    {
      subsymb.erase(po - nc,nc);
    }
  }else if(subsymb[0] == '.')
  {
    subsymb.erase(0,1);
  }

}

void Units::reversePowerSymbol(std::string &subsymbol)
{
  unsigned int curInter;
  curInter = subsymbol.find(".");
  if(subsymbol.find("/") < curInter)curInter = subsymbol.find("/");
  while(curInter < subsymbol.size())
  {
   (subsymbol[curInter] == '.')?
      subsymbol.replace(curInter,1,"/"):
      subsymbol.replace(curInter,1,".");
   curInter++;
   unsigned int tmp = subsymbol.find(".",curInter);
   if(subsymbol.find("/",curInter) < tmp)tmp = subsymbol.find("/",curInter);
   curInter = tmp;
  }
}

void Units::symbolToThePower(int r,const std::string &key)
{
  const std::string method("Units::symbolToThePower(int)");
  if(nDimensionOfUnits() == 0)return;
  std::string curUnit("");
  std::string tmpSymbol = contractedSymbol();
  for(unsigned int i = 0; i < tmpSymbol.size(); i++)
  {
    if(tmpSymbol[i] != '.' && tmpSymbol[i] != '/')curUnit += tmpSymbol[i];
    if(tmpSymbol[i] == '.' || tmpSymbol[i] == '/' || i == tmpSymbol.size() - 1 )
    {
      int nc(0);
      std::ostringstream po;
      int resultPower = getIntegerPower(parsePower(curUnit,nc),r,key); //power
      if(resultPower == 0)
      {
        symbol = "failed";
        return;
      }
      po << resultPower;
      std::string postr = po.str();
      if(resultPower == 1)postr = "";
      if(nc != 0)
      {
        (i != tmpSymbol.size() - 1)?tmpSymbol.replace(i - nc,nc,postr):
                                    tmpSymbol.replace(i - nc + 1,nc,postr); //resulting unit
      }else
      {
        (i != tmpSymbol.size() - 1)?tmpSymbol.insert(i,postr):tmpSymbol.insert(i + 1,postr);
      }
      curUnit.clear();
      i += postr.size();
    }
  }

  symbol = tmpSymbol;
}

int Units::getIntegerPower(int unit,int r, const std::string &key)
{
  if(key == "multiplication")
  {
     return unit * r;
  }else if(key == "division")
  {
    if(unit%r != 0)return 0;
    return unit / r;
  }else
  {
    const std::string method("int Units::getIntegerPower(int ,int , const std::string )");
    antiochError(method,"Key is not acceptable. This is a private method, there is a big problem...");
    return 0;
  }


}

double Units::FactorToSomeUnit(const Units &target) const
{
  if(isHomogeneous(target))
  {
    return getSIFactor()/target.getSIFactor();
  }else
  {
    const std::string method("double Units::FactorToSomeUnit(const Units &)const");
    antiochError(method,"Units are not homogeneous:\n\"" + symbol + "\" and \"" + target.getSymbol() + "\".");
    return -1.;
  }
}

double Units::TranslatorToSomeUnit(const Units & target)  const
{
  if(isHomogeneous(target))
  {
    return ((getSITranslator() - target.getSITranslator())/target.getSIFactor());
  }
  else
  {
    const std::string method("double Units::TranslatorToSomeUnit(const Units &)const");
    antiochError(method,"Units are not homogeneous.");
    return -1.;
  }
}

Converter Units::raise(const Converter &tbm,int power) const
{
  return Converter(std::pow(tbm.geta(),power),(power != 1)?0.:tbm.getb());
}

Converter Units::raise(const Converter &tbm,double power) const
{
  return Converter(std::pow(tbm.geta(),power),(power != 1.)?0.:tbm.getb());
}

const bool Units::isInSymb(char c) const
{

  return (c != '/' && 
          c != '.' && 
          c != '-' &&
          !isNumber(c));
}

void Units::showAll(std::ostream &out)
{
  out << "Unit description:" << std::endl;
  out << "name: "   << name << std::endl;
  out << "symbol: " << symbol << std::endl;
  out << "SI decomposition: " << power << std::endl;
  out << "SI converter: " << toSI << std::endl << std::endl;
}

Units& Units::operator=(const Units & rhs)
{
  if(this == &rhs){return *this;}
  name = rhs.getName();
  symbol = rhs.getSymbol();
  toSI.clear();
  power.clear();
  fillInPower(true);
  return *this;
}

Units & Units::operator+=(const Units & rhs)
{

// the name
  if(!rhs.getName().empty())name  += " " + rhs.getName();

// the symbol
  if(rhs.getSymbol().empty())return *this;
  if(!symbol.empty())
  {
    symbol  += ".(" + rhs.getSymbol() + ")";
  }else
  {
    symbol = rhs.getSymbol();
  }
  toSI  *= rhs.getSIcoef();
  power += rhs.getPower();
  return *this;
}

Units & Units::operator-=(const Units & rhs)
{
  if(!rhs.getName().empty())name += " / " + rhs.getName();
  if(!rhs.getSymbol().empty())
  {
    symbol  += "/(" + rhs.getSymbol() + ")";
  }
  developSymbol(symbol);
  toSI /= rhs.getSIcoef();
  power -= rhs.getPower();

  return *this;
}

Units & Units::operator*=(int r)
{
   power *= r;
   toSI = raise(toSI,r);
   symbolToThePower(r,"multiplication");
   if(symbol == "failed")
   {
     symbol = getSISymb();
   }

   return *this;
}

Units & Units::operator/=(int r)
{
  power /= r;//check consistency of root
  symbolToThePower(r,"division");
  if(symbol == "failed")
  {
     symbol = getSISymb();
  }
  toSI = raise(toSI,1./(double)r);

  return *this;
}


Units Units::operator+(const Units & rhs) const
{
  return (Units(*this) += rhs);
}

Units Units::operator-(const Units & rhs) const
{
  return (Units(*this) -= rhs);
}

Units Units::operator*(int r) const
{
  return (Units(*this) *= r);
}

Units Units::operator/(int r) const
{
  return (Units(*this) /= r);
}


// Public functions
std::string allKnownUnits()
{
  std::ostringstream out;
  out << std::left;
  out << std::setw(15)  << "Unit symbol";
  out << std::setw(25) << "Unit name";
  out << std::setw(20) << "Coefficient";   
  out << std::setw(50) << "Projection" << std::endl << "--" << std::endl;
  for(int i = 0; i < nKnownUnits; i++)
  {
     std::ostringstream symb,name,coef,power;
     symb  << knownUnits[i].getSymbol();
     name  << knownUnits[i].getName();
     coef  << knownUnits[i].getSIcoef();
     power << knownUnits[i].getPower();

     out << std::setw(15) << symb.str();
     out << std::setw(25) << name.str();
     out << std::setw(20) << coef.str();
     out << std::setw(50) << power.str() << std::endl;
  }

  return out.str();
}

std::string allKnownHomogeneousUnits(std::string target)
{
  std::ostringstream out;
  out << std::left;
  out << "Homogeneous known units to " << target << std::endl;
  out << std::setw(15)  << "Unit symbol";
  out << std::setw(25) << "Unit name";
  out << std::setw(20) << "Coefficient";   
  out << std::setw(50) << "Projection" << std::endl << "--" << std::endl;
  for(int i = 0; i < nKnownUnits; i++)
  {
     if(!knownUnits[i].isHomogeneous(target))continue;
     std::ostringstream symb,name,coef,power;
     symb  << knownUnits[i].getSymbol();
     name  << knownUnits[i].getName();
     coef  << knownUnits[i].getSIcoef();
     power << knownUnits[i].getPower();

     out << std::setw(15) << symb.str();
     out << std::setw(25) << name.str();
     out << std::setw(20) << coef.str();
     out << std::setw(50) << power.str() << std::endl;
  }

  return out.str();
}


std::string allKnownPrefixes()
{
  std::ostringstream out;
  out << std::left;
  out << std::setw(18)  << "Prefixe symbol";
  out << std::setw(18) << "Prefixe name";
  out << std::setw(15) << "Value" << std::endl << "--" << std::endl;
  for(int i = 0; i < nSIPrefixes; i++)
  {
     std::ostringstream symb,name,coef;
     symb  << Prefixes[i].symbol();
     name  << Prefixes[i].name();
     coef  << Prefixes[i].value();

     out << std::setw(18) << symb.str();
     out << std::setw(18) << name.str();
     out << std::setw(15) << coef.str() << std::endl;
  }

  return out.str();
}

void AreDefsOk()
{
  const std::string method("void AreDefsOk()");

  for(int i = 0; i < nKnownUnits; i++)
  {
     for(int p = 0; p < nSIPrefixes; p++)
     {
        std::string currentU = Prefixes[p].symbol() + knownUnits[i].getSymbol();
        for(int j = 0; j < nKnownUnits; j++)
        {
           for(int q = 0; q < nSIPrefixes; q++)
           {
             if(j == i && p == q)continue;
              std::string testU = Prefixes[q].symbol() + knownUnits[j].getSymbol();
              if(currentU == testU)
                antiochError(method,
                                "The file unit_defs.hpp is not good, check your knownUnits[] and Prefixes[].");
           }
        }
     }  
  }
}

}
