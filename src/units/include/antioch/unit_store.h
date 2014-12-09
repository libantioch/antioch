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
#ifndef ANTIOCH_UNITS_STORE_H
#define ANTIOCH_UNITS_STORE_H

/*!\file unit_store.hpp
 * \brief Provided units and prefixes
 *
 * This file contains the known units to be
 * used for unit management. Some physical constant
 * (double form) are also defined. The known units
 * begins with the SI basis:
 * \f[
 *      [\mathrm{meter,kilogram,second,ampere,kelvin,mole,candela,radian}]
 * \f]
 * 
 * The definition of a unit uses the fully-descriptive constructor
 * Units::UnitBase<T>(string sym,string na,double conva,double convb,int mi,int kgi=0, int si=0, int Ai=0, int Ki=0, int moli=0, int cdi=0, int radi=0);
 * Any unit is thus characterized as follow:
 *
 * UnitBase<T>("mmHg","millimeter of mercury",133.322387415,0.,-1,1,-2).
 *
 * The prefixes are then defined, in the form
 *
 * SIPrefixes<T>("mu",1e-6).
 *
 */

//Antioch
#include "antioch/unit_base.h"
#include "antioch/siprefix.h"
#include "antioch/math_constants.h"

//C++
#include <map>
#include <vector>


namespace Antioch{

namespace UnitBaseConstant{

  template <typename T = double>
  class UnitBaseStore{
        public:
    ~UnitBaseStore(){}
    UnitBaseStore()
    {
    /* base UnitBase, 
     * UnitBase<T>(string symbol, string name, 
     * double a, double b, 
     * int m, int kg = 0, int s = 0, int A = 0, int K = 0, int mol = 0, int cd = 0, int rad = 0)
     */
    store.push_back(UnitBase<T>("m",  "meter",   1.0L,0.,1));
    store.push_back(UnitBase<T>("kg", "kilogram",1.0L,0.,0,1));
    store.push_back(UnitBase<T>("s",  "second",  1.0L,0.,0,0,1));
    store.push_back(UnitBase<T>("A",  "ampere",  1.0L,0.,0,0,0,1));
    store.push_back(UnitBase<T>("K",  "kelvin",  1.0L,0.,0,0,0,0,1));
    store.push_back(UnitBase<T>("mol","mole",    1.0L,0.,0,0,0,0,0,1));
    store.push_back(UnitBase<T>("cd", "candela", 1.0L,0.,0,0,0,0,0,0,1));
    store.push_back(UnitBase<T>("rad","radian",  1.0L,0.,0,0,0,0,0,0,0,1));
    /* others known Units using constructor 
     * UnitBase<T>(string symbol, string name, 
     * double a, double b, 
     * int m, int kg=0, int s=0, int A=0, int K=0, int mol=0, int cd=0, int rad=0)
     */
    //length
    store.push_back(UnitBase<T>("in", "inch",             0.0254L,          0.,1));
    store.push_back(UnitBase<T>("ft", "foot",             0.3048L,          0.,1));
    store.push_back(UnitBase<T>("ua", "astronomical unit",1.49597870700e11L,0.,1)); //~ Sun-Earth distance, BIPM symbol
    store.push_back(UnitBase<T>("ang","angstrom",         1e-10L,           0.,1)); 
    //mass
    store.push_back(UnitBase<T>("Da","dalton",             1.660538921e-27L,0.,0,1));
    store.push_back(UnitBase<T>("u", "unified atomic mass",1.660538921e-27L,0.,0,1));
    //Time
    store.push_back(UnitBase<T>("min", "minute",60.0L,  0.,0,0,1));
    store.push_back(UnitBase<T>("hour","hour",  3600.0L,0.,0,0,1));
    //Current
    //Temperature
    store.push_back(UnitBase<T>("degC","degree Celsius",  1.L,    273.15L,          0,0,0,0,1));
    store.push_back(UnitBase<T>("degF","degree Farenheit",5.L/9.L,459.57L*5.0L/9.0L,0,0,0,0,1));
    //angle
    store.push_back(UnitBase<T>("deg", "degree",    Constants::pi<T>()/180.0L,              0.,0,0,0,0,0,0,0,1));
    store.push_back(UnitBase<T>("\'",  "Arcminute", Constants::pi<T>()/(180.0L*60.0L),      0.,0,0,0,0,0,0,0,1));
    store.push_back(UnitBase<T>("\'\'","Arcsecond", Constants::pi<T>()/(180.0L*60.0L*60.0L),0.,0,0,0,0,0,0,0,1));
    //volume
    store.push_back(UnitBase<T>("l","litre",1e-3L,0.,3));
    store.push_back(UnitBase<T>("L","litre",1e-3L,0.,3));
    //Force m.kg.s-2
    store.push_back(UnitBase<T>("N",  "newton",1.0L, 0.,1,1,-2));
    store.push_back(UnitBase<T>("dyn","dyne",  1e-5L,0.,1,1,-2));
    //Pressure N.m-2 (m-1.kg.s-2)
    store.push_back(UnitBase<T>("Pa",  "pascal",               1.0L,            0.,-1,1,-2));
    store.push_back(UnitBase<T>("bar", "bar",                  1e5L,            0.,-1,1,-2));
    store.push_back(UnitBase<T>("at",  "technical atmosphere", 9.80665e4L,      0.,-1,1,-2));
    store.push_back(UnitBase<T>("atm", "atmosphere",           1.01325e5L,      0.,-1,1,-2));
    store.push_back(UnitBase<T>("Torr","Torr",                 101325.0L/760.0L,0.,-1,1,-2)); // 1 atm / 760 by def
    store.push_back(UnitBase<T>("psi", "pound per square inch",6.895e3L,        0.,-1,1,-2));
    store.push_back(UnitBase<T>("mmHg","millimeter of mercury",133.322387415L,  0.,-1,1,-2));
    //viscosity Pa.s (m-1.kg.s-1)
    store.push_back(UnitBase<T>("P","poise",0.1L,0.,-1,1,-1));
    //Power W (m2.kg.s-3)
    store.push_back(UnitBase<T>("W","watt" ,1.0L,0.,2 ,1,-3));
    //Energy N.m (m2.kg.s-2)
    store.push_back(UnitBase<T>("J",    "joule",       1.0L,            0.,2,1,-2));
    store.push_back(UnitBase<T>("cal",  "calorie",     4.1868L,         0.,2,1,-2));//International table
    store.push_back(UnitBase<T>("calth","calorie thermodynamic", 4.184L,0.,2,1,-2));//thermodynamics
    store.push_back(UnitBase<T>("eV",   "electronVolt",1.602176565e-19L,0.,2,1,-2));
    store.push_back(UnitBase<T>("erg",  "erg",         1.e-7L,          0.,2,1,-2));
    store.push_back(UnitBase<T>("Ha",   "hartree",     4.35974434e-18L, 0.,2,1,-2));
    //Electric charge A.s
    store.push_back(UnitBase<T>("C","coulomb",1.0L,0.,0,0,1,1));
    //Frequency s-1
    store.push_back(UnitBase<T>("Hz","herzt",    1.0L,   0.,0,0,-1));
    store.push_back(UnitBase<T>("Ci","curie",    3.7e10L,0.,0,0,-1));
    store.push_back(UnitBase<T>("Bq","becquerel",1.0L,   0.,0,0,-1));
    // dipole moment C.m = A.s.m
    store.push_back(UnitBase<T>("D","debye",3.335641e-30L,0.,1,0,1,1)); // http://cccbdb.nist.gov/debye.asp)
    //no unit
    store.push_back(UnitBase<T>("molecule","molecule",1.0L,0.,0,0,0,0,0,0));
    store.push_back(UnitBase<T>("photon",  "photon",  1.0L,0.,0,0,0,0,0,0));
    
    _n_known_units = store.size();
    
      for(int i = 0; i < (int)store.size(); i++)
      {
        map_store[store[i].symbol()] = i;
      }
    }
    
    int stored_index(const std::string &symb)  const 
    {
      return (symb == "g")?map_store.at("kg"): //special case for gram
                                              (map_store.count(symb))?map_store.at(symb):-1;
    }
    
    const UnitBase<T> stored(const int &iunit) const
    {
      return store[iunit];
    }
    
    int n_known_units() const 
    {
      return _n_known_units;
    }

        private:
    std::map<std::string,unsigned int> map_store;
    std::vector<UnitBase<T> > store;
    unsigned int _n_known_units;
  };


/*!Prefixes, SI
 * micro is mu
 *
 *  - y,  1e-24: yocto
 *  - z,  1e-21: zepto
 *  - a,  1e-18: atto
 *  - f,  1e-15: femto
 *  - p,  1e-12: pico
 *  - n,  1e-9:  nano
 *  - mu, 1e-6:  micro
 *  - m,  1e-3:  milli
 *  - c,  1e-2:  centi
 *  - d,  1e-1:   deci
 *  - da, 1e1:   deca
 *  - h,  1e2:   hecto
 *  - k,  1e3:   kilo
 *  - M,  1e6:   mega
 *  - G,  1e9:   giga
 *  - T,  1e12:  tera
 *  - P,  1e15:  peta
 *  - E,  1e18:  exa
 *  - Z,  1e21:  zetta
 *  - Y,  1e24:  yotta
 */
template <typename T = double>
class SIPrefixeStore{
   public:
     ~SIPrefixeStore(){}
     SIPrefixeStore(){
     store.push_back(SIPrefixes<T>("y", "yocto",1e-24L));
     store.push_back(SIPrefixes<T>("z", "zepto",1e-21L));
     store.push_back(SIPrefixes<T>("a", "atto", 1e-18L));
     store.push_back(SIPrefixes<T>("f", "femto",1e-15L));
     store.push_back(SIPrefixes<T>("p", "pico", 1e-12L));
     store.push_back(SIPrefixes<T>("n", "nano", 1e-9L));
     store.push_back(SIPrefixes<T>("mu","micro",1e-6L));
     store.push_back(SIPrefixes<T>("m", "milli",1e-3L));
     store.push_back(SIPrefixes<T>("c", "centi",1e-2L));
     store.push_back(SIPrefixes<T>("d", "deci", 1e-1L));
     store.push_back(SIPrefixes<T>("da","deca", 1e1L));
     store.push_back(SIPrefixes<T>("h", "hecto",1e2L));
     store.push_back(SIPrefixes<T>("k", "kilo", 1e3L));
     store.push_back(SIPrefixes<T>("M", "mega", 1e6L));
     store.push_back(SIPrefixes<T>("G", "giga", 1e9L));
     store.push_back(SIPrefixes<T>("T", "tera", 1e12L));
     store.push_back(SIPrefixes<T>("P", "peta", 1e15L));
     store.push_back(SIPrefixes<T>("E", "exa",  1e18L));
     store.push_back(SIPrefixes<T>("Z", "zetta",1e21L));
     store.push_back(SIPrefixes<T>("Y", "yotta",1e24L));
     
     _n_prefixes = store.size();
     
       for(int i = 0; i < (int)store.size(); i++)
       {
           map_store[store[i].symbol()] = i;
       }
     }
     
     int stored_index(const std::string &symb) const 
     {
        return (map_store.count(symb))?map_store.at(symb):-1;
     }

     const SIPrefixes<T> stored(const int &ipre) const
     {
        return store[ipre];
     }

     int n_known_prefixes() const 
     {
        return _n_prefixes;
     }
  
     private:
       std::map<std::string, int> map_store;
       std::vector<SIPrefixes<T> > store;
       unsigned int _n_prefixes;
  };


}//end namespace UnitBase


}//end namespace Antioch
#endif
