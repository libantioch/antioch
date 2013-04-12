//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
#ifndef _UNITS_DEFINITIONS_
#define _UNITS_DEFINITIONS_

#include "antioch/Units.hpp"

/*!\file unit_defs.hpp
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
 * Units::Units(string sym,string na,double conva,double convb,int mi,int kgi=0, int si=0, int Ai=0, int Ki=0, int moli=0, int cdi=0, int radi=0);
 * Any unit is thus characterized as follow:
 *
 * Units("mmHg","millimeter of mercury",133.322387415,0.,-1,1,-2).
 *
 * The prefixes are then defined, in the form
 *
 * SIPrefixes("mu",1e-6).
 *
 */

namespace Antioch{

//some useful physical constants
namespace unitConstant{//
const double Av = 6.02214129e23; //Avogadro number to precision
const double pi = 3.141592653589793; //pi to 10^-15
}

const Units knownUnits[]={//
/* base Units, 
 * Units::Units(string symbol, string name, 
 * double converter a, double converter b, 
 * int m, int kg = 0, int s = 0, int A = 0, int K = 0, int mol = 0, int cd = 0, int rad = 0)
 */
Units("m","meter",1.,0.,1),
Units("kg","kilogram",1.,0.,0,1),
Units("s","second",1.,0.,0,0,1),
Units("A","ampere",1.,0.,0,0,0,1),
Units("K","kelvin",1.,0.,0,0,0,0,1),
Units("mol","mole",1.,0.,0,0,0,0,0,1),
Units("cd","candela",1.,0.,0,0,0,0,0,0,1),
Units("rad","radian",1.,0.,0,0,0,0,0,0,0,1),
/* others known Units using constructor 
 * Units::Units(string symbol, string name, 
 * double converter a, double converter b, 
 * int m, int kg=0, int s=0, int A=0, int K=0, int mol=0, int cd=0, int rad=0)
 */
//length
Units("in","inche",0.0254,0.,1),
Units("ft","foot",0.3048,0.,1),
Units("ua","astronomical unit",1.49597870700e11,0.,1), //~ Sun-Earth distance, BIPM symbol
Units("ang","angstrom",1e-10,0.,1), 
//mass
Units("Da","dalton",1.660538921e-27,0.,0,1),
Units("u","unified atomic mass",1.660538921e-27,0.,0,1),
//Time
Units("min","minute",60.,0.,0,0,1),
Units("hour","hour",3600.,0.,0,0,1),
//Current
//Temperature
Units("degC","degree Celsius",1.,273.15,0,0,0,0,1),
Units("degF","degree Farenheit",5./9.,459.57*5./9.,0,0,0,0,1),
//angle
Units("deg","degree",unitConstant::pi/180.,0.,0,0,0,0,0,0,0,1),
Units("\'","Arcminute",unitConstant::pi/(180.*60.),0.,0,0,0,0,0,0,0,1),
Units("\'\'","Arcsecond",unitConstant::pi/(180.*60.*60.),0.,0,0,0,0,0,0,0,1),
//volume
Units("l","litre",1e-3,0.,3),
Units("L","litre",1e-3,0.,3),
//Force m.kg.s-2
Units("N","newton",1.,0.,1,1,-2),
Units("dyn","dyne",1e-5,0.,1,1,-2),
//Pressure N.m-2 (m-1.kg.s-2)
Units("Pa","pascal",1.,0.,-1,1,-2),
Units("bar","bar",1e5,0.,-1,1,-2),
Units("at","technical atmosphere",9.80665e4,0.,-1,1,-2),
Units("atm","atmosphere",1.01325e5,0.,-1,1,-2),
Units("Torr","Torr",101325./760.,0.,-1,1,-2), // 1 atm / 760 by def
Units("psi","pound per square inch",6.895e3,0.,-1,1,-2),
Units("mmHg","millimeter of mercury",133.322387415,0.,-1,1,-2),
//Power W (m2.kg.s-3)
Units("W","watt",1.,0.,2,1,-3),
//Energy N.m (m2.kg.s-2)
Units("J","joule",1.,0.,2,1,-2),
Units("cal","calorie",4.184,0.,2,1,-2),
Units("eV","electronVolt",1.602176565e-19,0.,2,1,-2),
Units("erg","erg",1.e-7,0.,2,1,-2),
Units("Ha","hartree",4.35974434e-18,0.,2,1,-2),
//Electric charge A.s
Units("C","coulomb",1.,0.,0,0,1,1),
//Frequency s-1
Units("Hz","herzt",1.,0.,0,0,-1),
Units("Ci","curie",3.7e10,0.,0,0,-1),
Units("Bq","becquerel",1.,0.,0,0,-1),
//no unit
Units("molecule","molecule",1.0,0.,0,0,0,0,0,0),
Units("photon","photon",1.0,0.,0,0,0,0,0,0)
};

const int nKnownUnits = sizeof(knownUnits)/sizeof(knownUnits[0]);

/*!Prefixes, quasi-SI, centi and deci added,
 * micro is mu
 *
 *  - a, 1e-18: atto
 *  - f, 1e-15: femto
 *  - p, 1e-12: pico
 *  - n, 1e-9: nano
 *  - mu, 1e-6: micro
 *  - m, 1e-3: milli
 *  - c, 1e-2: centi
 *  - d, 1e1: deci
 *  - k, 1e3: kilo
 *  - M, 1e6: mega
 *  - G, 1e9: giga
 *  - T, 1e12: tera
 *  - P, 1e15: peta
 *  - E, 1e18: exa
 */
const SIPrefixes Prefixes[]={
SIPrefixes("a","atto",1e-18),
SIPrefixes("f","femto",1e-15),
SIPrefixes("p","pico",1e-12),
SIPrefixes("n","nano",1e-9),
SIPrefixes("mu","micro",1e-6),
SIPrefixes("m","milli",1e-3),
SIPrefixes("c","centi",1e-2),
SIPrefixes("d","deci",1e-1),
SIPrefixes("k","kilo",1e3),
SIPrefixes("M","mega",1e6),
SIPrefixes("G","giga",1e9),
SIPrefixes("T","tera",1e12),
SIPrefixes("P","peta",1e15),
SIPrefixes("E","exa",1e18)
};

const int nSIPrefixes = sizeof(Prefixes)/sizeof(Prefixes[0]);
}
#endif
