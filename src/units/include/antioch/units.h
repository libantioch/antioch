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
#ifndef ANTIOCH_UNITS_H
#define ANTIOCH_UNITS_H

//Antioch
#include "antioch/metaprogramming.h"
#include "antioch/antioch_asserts.h"
#include "antioch/siprefix.h"
#include "antioch/insi.h"
#include "antioch/converter.h"
#include "antioch/unit_defs.h"

//C++
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>

namespace Antioch{

/*!\file units.h
 * \brief Advanced unit class
 *
 * \class Units
 * \brief An advanced unit class
 *
 * This class contains:
 *   - a std::string for the symbol,
 *   - a std::string for the name,
 *   - a InSI object for the power vector
 *   - a Converter object for the coefficient
 *
 * The idea of this class is to keep track of the unit without
 * performing any operation on the values. To characterize an unit,
 * the unit is projected against the SI basis of units, with the
 * addition of angle unit:
 *
 * (meter,kilogram,second,ampere,kelvin,mol,candela,radian)
 *
 * Thus any unit can be expressed as a linear composition of
 * this std::vector. This is the role of the class InSI to store this
 * linear composition. This linear composition is called the power vector.
 *
 * With the linear composition, comes a coefficient that correspond
 * to the conversion of the considered unit to the SI system. The class
 * Converter takes care of storing it.
 *
 * Thus this class performs a parsing of the symbol of the given unit,
 * compares each unit of the symbol to a set of known units, and construct
 * the linear composition and coefficient by successive additions of
 * the units of the symbol. For instance, the unit m.A/s is seen as
 * \f$\mathrm{m} + \mathrm{A} - \mathrm{s}\f$, and thus the composition 
 * and coefficient are calculated correspondingly:
 *   - adding two units means adding the compositions and multiplying the coefficients
 *   \f[ 
 *         \mathrm{Units}_1 + \mathrm{Units}_2 \Rightarrow 
 *                                           \left\{\begin{array}{l}
 *                                                     \mathrm{InSI}_1 + \mathrm{InSI}_2\\
 *                                                     \mathrm{Converter}_1 \times \mathrm{Converter}_2
 *                                                   \end{array}\right.
 *   \f]
 *   - substracting two units means substracting the compositions and dividing the coefficients
 *   \f[ 
 *         \mathrm{Units}_1 - \mathrm{Units}_2 \Rightarrow 
 *                                           \left\{\begin{array}{l}
 *                                                     \mathrm{InSI}_1 - \mathrm{InSI}_2\\
 *                                                     \frac{\mathrm{Converter}_1}{\mathrm{Converter}_2}
 *                                                   \end{array}\right.
 *   \f]
 *
 * At each unit parsing level, the unit is parsed between a prefix, a symbol and a power.
 * The symbol "cm-3" will be parsed as a prefix "c", a symbol "m" and a power "-3". To this
 * parsing, the symbol before the unit is considered. If there is none or it is a dot, nothing
 * to do, but if this is a slash, the power has to be multiplied by -1. Thus
 * the Converter object will have the corresponding SIPrefixes values (Prefixe[6] is this case,
 * which is ("c",1e-2)), and the InSI power vector will add the value "-3" to the corresponding
 * power (0 in this case) if the character before the unit is a dot or does not exist,
 * the value 3 if this character is a slash. The correspondancies between the index in the power vector and
 * the unit is done explicitely.
 *
 * A special treatment is reserved for the kg unit, as it contains a prefix.
 *
 * This method supports a wide variety of symbol writing, like the use of parenthesises. The
 * used symbol is a fully developped one. 
 */
template <typename T = double>
class Units{
      public:

/*!\brief Fully descriptive constructor, mainly for internal purposes*/
        Units(const std::string &sym, const std::string &na,
                const T &conva, const T &convb,
                int mi,int kgi=0, int si=0, int Ai=0, int Ki=0, int moli=0, int cdi=0, int radi=0);//constructors for known Units
/*!\brief User-friendly constructor, from a symbol and an optionnal name*/
        Units(const std::string &sym, const std::string &na="");  //constructors when automatic conversion
/*!\brief Constructor for a user imposed conversion*/
        Units(const std::string & sym, const Converter<T> &conv, const std::string &na="");//constructor for given conversion
/*!\brief Default constructor*/
        Units(){}

/*!\brief Default destructor*/
        ~Units(){}

/*! \brief Comparison operator
 *
 * Compares the power vector and the coefficient, regardless
 * of the symbol and the name.
 */
        bool const operator== (const Units<T> &rhs) const;
/*! \brief Comparison operator, not equal is not "equal"*/
        bool const operator!= (const Units<T> &rhs) const {return !(*this == rhs);}
/*!\brief Assignement operator
 *
 * It assigns the symbol std::string,
 * the power vector and coefficient are filled.
 */
        Units<T> &operator= (const Units<T> & rhs);
/*! \brief Adding operator
 *
 * Adding units occurs when multiplying physical quantities.
 *
 * Adding a unit needs to add first the name if the added
 * name is not empty. Addition of name consist simply of
 * putting them side by side with a blank space in between.
 * 
 * The symbol is added in the same fashion, but with a dot
 * in between to ensure parsability. If the current symbol is
 * empty, them it is simply copied.
 * 
 * The power is added and the
 * coefficient is multiplied.
 */
        Units<T> &operator+=(const Units<T> & rhs);
/*!\brief Substracting operator
 *
 * Substracting units occurs when dividing physical quantities.
 *
 * Substracting is close to adding (see Units &operator+=(const Units&)),
 * with little differences. The adding of the name is done
 * with the sequence blank-slash-blank (" / ") between the names.
 * The symbol needs a bit more processing: if the added symbol is
 * not empty, is it put between parenthesises and added with a slash
 * in between. Then the void develop_symbol() method is called to
 * suppress the parenthesises and obtain a correct symbol.
 *
 * The coefficient is divided and
 * the power vector substracted.
 */
        Units<T> &operator-=(const Units<T> & rhs);
/*!\brief Adding operator, uses Units &operator+=(const Units&)*/
        Units<T>  operator+ (const Units<T> & rhs) const;
/*!\brief Substracting operator, uses Units &operator-=(const Units&)*/
        Units<T>  operator- (const Units<T> & rhs) const;
/*!\brief Multiplying operator by an integer
 *
 * This operator is needed for power elevation operation of a physical quantity.
 *
 * The operations needed are as follow:
 *   - parsing the symbol so all the power are multiplied by the integer,
 *   done by the method void symbol_to_the_power(int,std::string);
 *   - multiplying the power vector by the integer;
 *   - elevating the coefficient to the power, done by the method
 *   Converter raise(const Converter&, int).
 *
 * Elevating to the power with a double is not supported, and the physical sense
 * of it anyway is not obvious...
 */
        Units<T> &operator*=(int r);
/*!\brief Multiplying operator by an integer, uses Units &operator*=(int)*/
        Units<T>  operator* (int r)             const;
/*!\brief dividing operator by an integer
 *
 * This operator is needed for root operation of a physical quantity.
 *
 * The operations needed are as follow:
 *   - dividing the power vector by the integer, test that the resulting
 *     power is an integer;
 *   - parsing the symbol so all the power are divided by the integer,
 *     done by the method void symbol_to_the_power(int,std::string),
 *     if the symbol contains unit symbol that gives non integer
 *     resulting power, then the SI symbol is called;
 *   - dividing the coefficient to the power, done by the method
 *     Converter raise(const Converter&, double).
 */
        Units<T> &operator/=(int r);
/*!\brief Dividing operator by an integer, uses Units &operator/=(int)*/
        Units<T>  operator/ (int r)             const;

/*!\brief Alternative call to Units& operator=(const Units &) */
        void equalize(const Units<T> &rhs)      {*this  = rhs;}
/*!\brief Alternative call to Units& operator=(const Units &) */
        void equalize(Units<T> *rhs)            {*this  = *rhs;}
/*!\brief Alternative call to Units& operator=(const Units &) 
 *
 * This method should not really be used. It is provided for simplicity,
 * but equalizing a Units object to a std::string is weird. Instead, the 
 * user should prefer to write ``*this = Units(std::string)''. This is exactly
 * what does the method.
 *
 * It is provided to ensure that the idea of equalizing the unit to
 * a symbol is coded, but in this object, a unit is more than just
 * a symbol.
 */
        void equalize(std::string unit)       {*this  = Units<T>(unit);}
/*!\brief Alternative call to Units& operator+=(const Units &) */
        void add(const Units<T> &rhs)            {*this += rhs;}
/*!\brief Alternative call to Units& operator+=(const Units &) */
        void add(Units<T> *rhs)                  {*this += *rhs;}
/*!\brief Alternative call to Units& operator+=(const Units &) 
 *
 * Same remarks than void equalize(std::string).
 */
        void add(std::string unit)            {*this += Units<T>(unit);}
/*!\brief Alternative call to Units& operator-=(const Units &)*/
        void substract(const Units<T> &rhs)      {*this -= rhs;}
/*!\brief Alternative call to Units& operator-=(const Units &)*/
        void substract(Units<T> *rhs)            {*this -= *rhs;}
/*!\brief Alternative call to Units& operator-=(const Units &) 
 *
 * Same remarks than void equalize(std::string).
 */
        void substract(const std::string &unit)      {*this -= Units<T>(unit);}
/*!\brief Alternative method to Units& operator/=(int)*/
        void root(int r)                      {*this /= r;}

/*! \brief Units setter, sets the unit to the given symbol and name*/
        void set_unit(const std::string &sym, std::string na);
/*! \brief Units setter, sets the unit to the given symbol and name*/
        void set_unit(const std::string &sym)       {set_unit(sym,"");}
/*! \brief Name setter*/
        void set_name(const std::string &na)        {name = na;}
/*! \brief Symbol setter*/
        void set_symbol(const std::string &symb)    {symbol = symb;}
/*! \brief Coefficient setter*/
        void set_SI_converter(const Converter<T> &conv) {Antioch::init_clone(toSI,conv);}
/*! \brief Power std::vector setter*/
        void set_power_to_SI(const InSI &mat)         {power = mat;}

/*!\brief Homogenity testing with another Units.
 *
 * Testing the homogeneity between units is different from testing
 * equality. J and cal are homogeneous but not equal for instance.
 * Thus this method merely compares the power vectors, if they're equal
 * the units are homogeneous, if not, they're not.
 */
        const bool is_homogeneous(const Units<T> &rhs) const;
/*!\brief Homogenity testing with a std::string.
 *
 * Testing homogeneity with a std::string needs to build
 * a Units object and testing homogeneity with it. In the
 * case of an empty std::string, instead of building a Units object
 * the boolean is_united() is tested.
 */
        const bool is_homogeneous(std::string target) const;
/*!\brief Test if the unit is non empty*/
        const bool is_united() const {return !power.empty();}


/*!\brief Calculates the factor to any given unit.
 *
 * Before calculating the factor to another unit, the
 * homogeneity is tested. If the units are homogeneous,
 * then the factor is the current factor divided by the target
 * factor:
 * \f[
 *    \left.\begin{array}{l}
 *              \mathrm{SI} = f_1 \times \mathrm{unit}_1 + t_1 \\
 *              \mathrm{SI} = f_2 \times \mathrm{unit}_2 + t_2
 *          \end{array}\right\}
 *    \Rightarrow \mathrm{unit}_2 = \mathrm{unit}_1 \times \frac{f_1}{f_2} + \frac{t_1 - t_2}{f_2}
 * \f]
 * It is important to understand that this is the factor \e to the wanted unit. Thus to
 * express a value of unit\f$_1\f$ in unit\f$_2\f$, one should use the conversion rule
 * given above. factor_to_some_unit(unit\f$_2\f$) will provide the factor \f$\frac{f_1}{f_2}\f$
 * while translator_to_some_unit(unit\f$_2\f$) will provide the translationnal term \f$\frac{t_1 - t_2}{f_2}\f$.
 */
        T factor_to_some_unit(const Units<T> & target)     const;
/*!\brief Calculates the factor to any given unit, uses factor_to_some_unit(const Units&).*/
        T factor_to_some_unit(const std::string &target)    const;
/*!\brief Calculates the translator to any given unit, see factor_to_some_unit(const Units&) for the equations.*/
        T translator_to_some_unit(const Units<T> & target)  const;
/*!\brief Calculates the translator to any given unit, uses translator_to_some_unit(const Units&).*/
        T translator_to_some_unit(const std::string &target)  const {return translator_to_some_unit(Units<T>(target));}

/*! \brief String name getter*/
        const std::string get_name()            const {return name;}
/*! \brief String symbol getter*/
        const std::string get_symbol()          const {return symbol;}
/*! \brief Coefficient getter*/
        const Converter<T> get_SI_coef() const {return toSI;}
/*! \brief Multiplicative coefficient getter*/
        T get_SI_factor()                const {return toSI.geta();}
/*! \brief Translationnal coefficient getter*/
        T get_SI_translator()            const {return toSI.getb();}
/*! \brief Power std::vector getter*/
        const InSI get_power()                  const {return power;}
/*! \brief Corresponding SI symbol getter
 *
 * This method returns the symbol that correspond to the power vector.
 * For instance a unit J/s will be returned as m2.kg.s-3
 */
        const std::string get_SI_symb()     const;
/*! \brief Power of the asked SI unit
 *
 * This method returns the power to the 
 * asked unit SI.
 */
        int get_SI_power(const std::string &SIask) const;
/*! \brief Corresponding SI symbol getter
 *
 * This method returns the symbol that of the SI unit (so a
 * coefficient of (1.,0.)). If none is found, it returns
 * const std::string get_SI_symb().
 * For instance a unit J/s will be returned as W.
 *
 * This method is not able to combine different SI units.
 * The unit J.s will not be recognized as J.s, only as
 * m2.kg.s-1
 */
        const std::string get_SI_convenient_symb() const;
/*! \brief Simplify the symbol when possible.
 *
 * This method will densify the symbol. As calculations
 * merely add the symbol, it is necessary to be able
 * to group together the same units and eliminate those
 * that compensate.
 *
 * The method will scan unit by unit the symbol std::string
 * and combine the units that are homogeneous. The unit "l"
 * is converted to "m3" to ensure meter harmonization. The chosen
 * unit symbol between two homogeneous units is the first symbol 
 * encountered.
 *
 * \todo The litre management is a fix specific to litre. What should
 * be done would be a pre-treatment of the symbol, test if its
 * power vector is a multiple of one of the unity of the base std::vector,
 * and make the change.
 */
        const std::string harmonized_symbol(const std::string &input = "")    const {return manipulate_symbol(input,false);}
/*! \brief Contract the symbol when possible.
 *
 * The contraction will aggregate only exactly equal
 * units.
 */
        const std::string contracted_symbol(const std::string &input = "")    const {return manipulate_symbol(input,true);}

/*! \brief Clear the unit.*/
        void clear(){symbol.clear();name.clear();toSI.clear();power.clear();}

/*! \brief showAll().*/
        void showAll(std::ostream &out = std::cout);

      private:

/*! \brief Root method for contracting or harmonizing the symbol
 *
 *  The symbol MUST be developped prior to be used in this method.
 */
        const std::string manipulate_symbol(std::string input, bool contract) const;


/*!\brief This method fills the power vector and will calculate the coefficient
 * if the bool doConv is set to true.
 *
 * The method scans the symbol, and calls bool parse_single_unit(int,std::string,bool)
 * for each unit found. The character in front of the unit ('.' or '/') is taken
 * into account into the first integer (resp. 1 or -1).
 * This is this latter method that actually calculates the coefficients and
 * power vector.
 */
        void fill_in_power(bool doConv);
/*!\brief Calculates the corresponding coefficient and power of the given unit.
 *
 * \param int signe, 1 or -1. Image of the sign before the unit '.' or '/'.
 * \param std::string unit. The unit to be parsed.
 * \param bool doConv. Switch to perform or not the coefficient calculation.
 *
 * The method works in two times:
 *   - first parsing the symbol into prefix + symbol + power
 *   - second update the power vector and, if doConv is true, the coefficient.
 * 
 * The parsing uses the methods int parse_power(std::string,int&) and parse_prefix_unit(int&,int&,std::string).
 * The idea is to scan the knownUnits[] and Prefixes[] array and find the corresponding unit
 * and prefix. The method int parse_power returns the power, while
 * the integer iPre and iUnit characterize the prefix and unit, with a value of -1 if
 * not found. An unknown unit will return an error.
 *
 * The update then is twofold. Adding to the power vector the power vector of
 * the found unit multiplied by the signe and the power:
 * \f[
       <p> = <p> + \left(<p>_{iUnit} \times s_u \times p_u\right)
 * \f]
 * with \f$<p>\f$ being the updated power vector, \f$p_u\f$ the power of the parsed unit,
 * \f$<p>_{iUnit}\f$ the power vector of the found knownUnit[],
 * $s_u$ the integer associated to the sign of the parsed unit.
 *
 * If the conversion is to be done, multiplying the coefficient by the coefficient 
 * of the known unit found, multiplied by the
 * coefficient of the prefix, raise to the power by the parsed power, updated by the signe:
 * \f[
 *      (c) = (c) \times \left( pr_{iPre} \times (c)_{iUnit}\right)^{s_u \times p_u}
 * \f]
 * with \f$(c)\f$ being the updated complex coefficient, \f$(c)_{iUnit}\f$ the complex
 * coefficient of the found unit, \f$pr_{iPre}\f$ the coefficient of the
 * found Prefixes (if found, if not, it is unity), 
 * \f$s_u\f$ the sign of the parsed unit and \f$p_u\f$ the
 * power of the parsed unit.
 */
        bool parse_single_unit(int signe, std::string unit, bool doConv);
/*!\brief Unit parser, parse a std::string of type Prefixes[] knownUnits[].
 *
 * \param int &iUnit. Output, indice of the found unit, if not found, set to -1.
 * \param int &iPre. Output, indice of the found prefixes, if not found, set to -1.
 * \param std::string unit. Input, std::string to be parsed.
 *
 * This method will parse a std::string, assuming the std::string is of the form
 * Prefixes[]knownUnits[]. No power should be present (see int parse_power(std::string,int&)
 * for power management).
 *
 * This method scans the Prefixes[] array, if one is found at the start of the
 * std::string, then the remaining std::string is tested to see if it is a unit (method
 * int indexUnit(std::string)). If the test fail the
 * scanning resumes. If a Prefixes[] and a knownUnits[] are found, then its over.
 * If no Prefixes[] associated to a knownUnits[] is found, then the whole std::string
 * is tested to be a knownUnits[].
 *
 * This method relies on the hypothesis that all possible pairs of Prefixes[] and knownUnits[]
 * (including no prefix) are different.
 */
        void parse_prefix_unit(int &iUnit,int &iPre, const std::string &unit)  const;
/*!\brief Scanner of knownUnits[], sends back the corresponding iUnit, -1 if not found.*/
        int indexUnit(std::string unit)                               const;
/*!brief parser for std::string of type "some_unit power_associated".
 *
 * \param std::string unit. Input, std::string to be parsed.
 * \param int &nc. Output, number of charaters of the power.
 * \return the power associated to the unit.
 */
        int parse_power(std::string unit, int &nc)                    const;
/*!\brief Small method to check if a character is a numerical.
 *
 * This method uses the values of the ascii table, which basically
 * means that is checks if the ascii integer associated to the
 * character is in the range [48,57].
 */
        bool is_number(char c)                                   const {return ((int)c <= 57 && (int)c >= 48);} //ascii table

/*!\brief Raise to an integer power a Converter object.
 *
 * The rising to the power of a type Converter is done as
 * follow:
 * \f[
 *     (a,b)^i =\left\{\begin{array}{ll} 
 *                      (a^i,b) & \mathrm{if}\ i\ \mathrm{equal\ 1} \\
 *                      (a^i,0) & \mathrm{else}
 *                      \end{array}\right.
 * \f]
 * \return A Converter object which is the results of the operation.
 */
        Converter<T> raise(const Converter<T> &tbm,int power)       const;
/*!\brief Overload to simplify division
 *
 * Another method, to secure multiplication. Division is thus
 * put apart, to secure against non integer power.
 */
        Converter<T> raise(const Converter<T> &tbm,double power)       const;
/*!\brief Supress the parenthesises
 *
 * Two types of parenthesises are to be treated:
 *   - serial parenthesises ().()
 *   - intricated parenthesises (())
 *
 * The idea is to treat every serial pairs one by one.
 * Thus every pair found is
 * sent back to develop_symbol(std::string&). Once the
 * we have treated the inside string by this recursive method,
 * we treat the string:
 *   - if the character before the opening parenthesis is '.',
 *      we just supress the parenthesises;
 *   - if it is '/', we reverse_power_symbol(std::string &) the string
 */
        void develop_symbol(std::string &subsymbol);
/*!\brief Treat part of symbol depending on parenthesises
 *
 * Inverse the power sign determinant ('.' or '/'). Note that
 * the first power will not be affected. Ex:
 * J.mol-1/K will be changed to J/mol-1.K, leaving
 * J unchanged.
 */
        void reverse_power_symbol(std::string &subsymbol);
/*!\brief Small method that checks if the symbol is within a unit std::string.
 *
 * This method checks if the symbol is not a number (bool is_number(char)),
 * and not any of the following character: '/', '.' and '-'. Then this
 * symbol is consider to belong to a unit symbol.
 */
        const bool is_in_symb(char c)               const;
/*!\brief Number of dimensions of the unit in the SI basis.*/
        const int n_dimension_of_units()             const;
/*!\brief Small method to add the SI symbol and power to a std::string.
 *
 * This method is useful in the std::string const Units::get_SI_symb() const
 * method.
 */
        const std::string add_SI(int pow,std::string SIsymb) const;
/*!\brief Rewrite the symbol to put it at a given power.
 *
 * This method scans the symbol, parse the units, and changes
 * the power by the power multiplied or divided by the given integer.
 * It uses the method get_integer_power(int,int,int) to ensure
 * integer power.
 */
        void symbol_to_the_power(int r, const int &key);
/*!\brief Multiply or divide a power
 *
 * It ensures that the result of the operation performed
 * is an integer. Natural in case of multiplication, to be
 * checked in case of division. It returns 0 if it is
 * not correct, this \e has to be checked by the wrapping
 * method.
 * The key is to choose between division or multiplication.
 * Key is:
 *   - 1 for multiplication
 *   - -1 for division
 */
        int get_integer_power(int unit,int r, const int &key);
/*!\brief Fill the name if empty and name obvious
 *
 * Usually we give only the symbol, even for simple unit.
 * Thus this method will fill the name if the name is obvious,
 * i.e. if the dimension of the power vector is 1. The
 * suffixe is given by the factor value.
 */
        void check_if_name(std::string &nameOut, const Units<T> &un) const;

/*!\brief Strings for the symbol.*/
        std::string symbol;
/*!\brief Strings for the name.*/
        std::string name;
/*!\brief A Converter for the coefficient.*/
        Converter<T> toSI;
/*!\brief An InSI for the power vector.*/
        InSI power;

};


template<typename T>
inline
Units<T>::Units(const std::string &sym, const std::string &na,
             const T &conva, const T &convb,
             int mi,int kgi, int si, int Ai, int Ki, int moli, int cdi, int radi): //fully descriptive constructor
  symbol(sym),
  name(na),
  toSI(conva,convb),
  power(mi,kgi,si,Ai,Ki,moli,cdi,radi)
{
  develop_symbol(symbol);
}

template<typename T>
inline
Units<T>::Units(const std::string &sym, const std::string &na):  //constructors when automatic conversion
       symbol(sym),
       name(na)
{
  develop_symbol(symbol);
  fill_in_power(true);
}

template<typename T>
inline
Units<T>::Units(const std::string &sym, const Converter<T> &conv, const std::string &na)://constructor for given conversion
       symbol(sym),
       name(na),
       toSI(conv)
{
  develop_symbol(symbol);
  fill_in_power(false);
}

template <typename T>
void Units<T>::fill_in_power(bool doConv)
{
  if(symbol.empty())return; // no unity
  std::string tmp(""),symboltmp(symbol);
  int signe(1),istart(0);

  while(symboltmp != contracted_symbol(symboltmp))symboltmp = contracted_symbol(symboltmp);
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
      if(!parse_single_unit(signe,tmp,doConv))
      {
        antioch_unit_error("In symbol " + symboltmp + ", unit \"" + tmp + "\" does not ring a bell");
        break;
      }
      signe = 1;
      tmp.clear();
      }else if(symboltmp[i] == '/')
      {
       if(!parse_single_unit(signe,tmp,doConv))
       {
         antioch_unit_error("In symbol " + symboltmp + ", unit \"" + tmp + "\" does not ring a bell");
         break;
       }
       signe = -1;
       tmp.clear();
       }else
       {
         tmp += symboltmp[i];
       }
   }
   if(!parse_single_unit(signe,tmp,doConv))
      antioch_unit_error("In symbol " + symboltmp + ", unit \"" + tmp + "\" does not ring a bell");
}

template <typename T>
void Units<T>::develop_symbol(std::string &subsymb)
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
        develop_symbol(insideStr);
        if(po != 0)if(subsymb[po - 1] == '/')reverse_power_symbol(insideStr);
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
    int nc(0),pow = - parse_power(curUnit,nc); //power
    if(pow != 1)
    {
      std::ostringstream np;
      np << pow;
      if(po < subsymb.size())
      {
         subsymb.replace(po - nc,nc,np.str());
      }else
      {
        subsymb += np.str();
      }
    }else
    {
      subsymb.erase(po - nc,nc);
    }
  }else if(subsymb[0] == '.')
  {
    subsymb.erase(0,1);
  }

}

template <typename T>
T Units<T>::factor_to_some_unit(const Units<T> &target) const
{
  if(is_homogeneous(target))
  {
    return get_SI_factor()/target.get_SI_factor();
  }else
  {
    antioch_unit_error("Units are not homogeneous:\n\"" + symbol + "\" and \"" + target.get_symbol() + "\".");
    return -1.;
  }
}

template<typename T>
T Units<T>::factor_to_some_unit(const std::string &target) const 
{
  return factor_to_some_unit(Units<T>(target));
}

template <typename T>
T Units<T>::translator_to_some_unit(const Units<T> & target)  const
{
  if(is_homogeneous(target))
  {
    return ((this->get_SI_translator() - target.get_SI_translator())/target.get_SI_factor());
  }
  else
  {
    antioch_unit_error("Units are not homogeneous.");
    return -1.;
  }
}

template <typename T>
void Units<T>::parse_prefix_unit(int &iUnit,int &iPre,const std::string& unit) const
{
  iPre = -1;
  iUnit = UnitBaseStorage::known_units().stored_index(unit);
//prefix, if exists, is one or two character
   if(iUnit == -1 && unit.size() > 1)
   {
     std::string pre = unit.substr(0,1);
     std::string un  = unit.substr(1,std::string::npos);
     iPre = UnitBaseStorage::known_prefixes().stored_index(pre);
     iUnit = UnitBaseStorage::known_units().stored_index(un);
     if(iPre == -1 && iUnit == -1 && unit.size() > 2)
     {
       pre = unit.substr(0,2);
       un  = unit.substr(2,std::string::npos);
       iPre = UnitBaseStorage::known_prefixes().stored_index(pre);
       iUnit = UnitBaseStorage::known_units().stored_index(un);
     }
   }
   return;
}

template <typename T>
bool Units<T>::parse_single_unit(int signe,std::string unit,bool doConv)
{
  int iUnit = -1, nc = 0;
  int ipower = parse_power(unit,nc); // find power
  int iPre = -1;

  unit = unit.substr(0,unit.length()-nc); //unit without power

  parse_prefix_unit(iUnit,iPre,unit);//find prefix and unit

  if(iUnit == -1)return false; //not found

  T pre = 1.;
  if(iPre != -1)pre = UnitBaseStorage::known_prefixes().stored(iPre).value<T>();
  if(UnitBaseStorage::known_units().stored(iUnit).symbol() == "kg")pre *= 1e-3;
  InSI powerTmp = UnitBaseStorage::known_units().stored(iUnit).power_array();
  power += (powerTmp * signe * ipower);
  if(doConv)
  {
    Converter<T> convTmp = UnitBaseStorage::known_units().stored(iUnit).converter() * pre;
    toSI += raise(convTmp,signe * ipower);
  }

  return true;
}

template <typename T>
int Units<T>::indexUnit(std::string unit) const
{
  int iUnit;

  if(unit.empty())return -1;

  if(unit == "g")//special case to adapt SI
  {
    unit = "kg";
  }
  for (iUnit = 0; iUnit < UnitBaseStorage::known_units().n_known_units() ; iUnit++)
  {
    if(unit == UnitBaseStorage::known_units().stored(iUnit).symbol())break;
  }

  if(iUnit >= UnitBaseStorage::known_units().n_known_units())iUnit=-1;


  return iUnit;
}

template <typename T>
int Units<T>::parse_power(std::string unit,int &nc) const
{
  int ip = 1, loc = unit.length();
  char c = unit[loc-1];
  nc = 0;
  while(is_number(c)){
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

  if(ip == 0)antioch_unit_error("Invalid power found: " + unit);

  return ip;
}

template <typename T>
const bool Units<T>::is_homogeneous(const Units<T> &rhs) const 
{
  return (power == rhs.get_power());
}

template <typename T>
const bool Units<T>::is_homogeneous(std::string target) const
{
  if(target.empty())
  {
    return (!this->is_united());
  }else
  {
    Units<T> rhs(target);
    return is_homogeneous(rhs);
  }
  
}

template <typename T>
int Units<T>::get_SI_power(const std::string &SIask) const
{
  if(SIask == "m")  return power.get_m();
  if(SIask == "kg") return power.get_kg();
  if(SIask == "s")  return power.get_s();
  if(SIask == "A")  return power.get_A();
  if(SIask == "K")  return power.get_K();
  if(SIask == "mol")return power.get_mol();
  if(SIask == "cd") return power.get_cd();
  if(SIask == "rad")return power.get_rad();

  antioch_unit_error(SIask + " is not a SI symbol");

  return 0;
}

template <typename T>
std::string const Units<T>::get_SI_symb() const
{
 std::string SISymb;
 SISymb.clear();

 if(power.get_m() != 0)SISymb += add_SI(power.get_m(),"m");

 if(power.get_kg() != 0)
 {
   if(!SISymb.empty())SISymb += ".";
   SISymb += add_SI(power.get_kg(),"kg");
 }
 if(power.get_s() != 0)
 {
   if(!SISymb.empty())SISymb += ".";
   SISymb += add_SI(power.get_s(),"s");
 }
 if(power.get_A() != 0)
 {
   if(!SISymb.empty())SISymb += ".";
   SISymb += add_SI(power.get_A(),"A");
 }
 if(power.get_K() != 0)
 {
   if(!SISymb.empty())SISymb += ".";
   SISymb += add_SI(power.get_K(),"K");
 }
 if(power.get_mol() != 0)
 {
   if(!SISymb.empty())SISymb += ".";
   SISymb += add_SI(power.get_mol(),"mol");
 }
 if(power.get_cd() != 0)
 {
   if(!SISymb.empty())SISymb += ".";
   SISymb += add_SI(power.get_cd(),"cd");
 }
 if(power.get_rad() != 0)
 {
   if(!SISymb.empty())SISymb += ".";
   SISymb += add_SI(power.get_rad(),"rad");
 }

 return SISymb;

}

template <typename T>
const std::string Units<T>::add_SI(int pow,std::string symb) const
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

template <typename T>
std::string const Units<T>::manipulate_symbol(std::string input, bool contract) const
{
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
      int curPower = parse_power(curUnit,nc); //power
//      strPower = curUnit.substr(curUnit.length() - nc,std::string::npos); //resulting power
      curUnit = curUnit.substr(0,curUnit.length() - nc); //resulting unit
      parse_prefix_unit(iUnit,iPre,curUnit);
 // unit is UnitBaseStorage::Prefixes[iPre].symbol() + UnitBaseStorage::knownUnits[iUnit].get_symbol() + strPower
      if(iUnit == -1)
      {
        antioch_unit_error("The unit \"" + curUnit + "\" is unknown. No harmonized symbol will be produced.");
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
         Units<T> tmp(unitvec[j]);
         bool same(tmp.get_symbol() == curUnit);//only contraction: strong condition
         if(!contract && !same)same = tmp.is_homogeneous(curUnit); //harmonizing: weak condition
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

template <typename T>
const int Units<T>::n_dimension_of_units() const
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

template <typename T>
const std::string Units<T>::get_SI_convenient_symb() const
{
  if(this->n_dimension_of_units() == 0)
  {
    std::string outStr("");
    return outStr;
  }else
  {
    for(int iUnit = 1; iUnit < UnitBaseStorage::known_units().n_known_units(); iUnit++)
    {
      if(is_homogeneous(UnitBaseStorage::known_units().stored(iUnit).symbol()) && 
         UnitBaseStorage::known_units().stored(iUnit).converter().geta() == 1. && 
         UnitBaseStorage::known_units().stored(iUnit).converter().getb() == 0.)return UnitBaseStorage::known_units().stored(iUnit).symbol();
    }
    return get_SI_symb();
  }
}

template <typename T>
void Units<T>::set_unit(const std::string &sym, std::string na)
{
  symbol = sym;
  develop_symbol(symbol);
  name = na;
  toSI.clear();
  power.clear();
  fill_in_power(true);
}

template <typename T>
void Units<T>::reverse_power_symbol(std::string &subsymbol)
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

template <typename T>
void Units<T>::symbol_to_the_power(int r,const int &key)
{
  if(n_dimension_of_units() == 0)return;
  std::string curUnit("");
  std::string tmpSymbol = contracted_symbol();
  for(unsigned int i = 0; i < tmpSymbol.size(); i++)
  {
    if(tmpSymbol[i] != '.' && tmpSymbol[i] != '/')curUnit += tmpSymbol[i];
    if(tmpSymbol[i] == '.' || tmpSymbol[i] == '/' || i == tmpSymbol.size() - 1 )
    {
      int nc(0);
      std::ostringstream po;
      int resultPower = get_integer_power(parse_power(curUnit,nc),r,key); //power
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

template <typename T>
int Units<T>::get_integer_power(int unit,int r, const int &key)
{
  if(key == 1)//multiplication
  {
     return unit * r;
  }else if(key == -1)//division
  {
    if(unit%r != 0)return 0;
    return unit / r;
  }else
  {
    std::cerr << "Key is not acceptable. This is a private method, there is a big problem..." << std::endl;
    antioch_error();
    return 0;
  }


}

template <typename T>
Converter<T> Units<T>::raise(const Converter<T> &tbm,int power) const
{
  return Converter<T>(std::pow(tbm.geta(),power),(power != 1)?0.:tbm.getb());
}

template <typename T>
Converter<T> Units<T>::raise(const Converter<T> &tbm, double power) const
{
  return Converter<T>(std::pow(tbm.geta(),power),(power != 1.)?0.:tbm.getb());
}

template <typename T>
const bool Units<T>::is_in_symb(char c) const
{

  return (c != '/' && 
          c != '.' && 
          c != '-' &&
          !this->is_number(c));
}

template <typename T>
void Units<T>::showAll(std::ostream &out)
{
  out << "Unit description:"  << std::endl;
  out << "name: "             << name   << std::endl;
  out << "symbol: "           << symbol << std::endl;
  out << "SI decomposition: " << power  << std::endl;
  out << "SI converter: "     << toSI   << std::endl << std::endl;
}

template <typename T>
Units<T>& Units<T>::operator=(const Units<T> & rhs)
{
  if(this == &rhs){return *this;}
  name = rhs.get_name();
  symbol = rhs.get_symbol();
  toSI.clear();
  power.clear();
  fill_in_power(true);
  return *this;
}

template <typename T>
Units<T> & Units<T>::operator+=(const Units<T> & rhs)
{

// the name
  if(!rhs.get_name().empty())name  += " " + rhs.get_name();

// the symbol
  if(rhs.get_symbol().empty())return *this;
  if(!symbol.empty())
  {
    symbol  += ".(" + rhs.get_symbol() + ")";
  }else
  {
    symbol = rhs.get_symbol();
  }
  toSI  *= rhs.get_SI_coef();
  power += rhs.get_power();
  return *this;
}

template <typename T>
Units<T> & Units<T>::operator-=(const Units<T> & rhs)
{
  if(!rhs.get_name().empty())name += " / " + rhs.get_name();
  if(!rhs.get_symbol().empty())
  {
    symbol  += "/(" + rhs.get_symbol() + ")";
  }
  this->develop_symbol(symbol);
  toSI  /= rhs.get_SI_coef();
  power -= rhs.get_power();

  return *this;
}

template <typename T>
Units<T> & Units<T>::operator*=(int r)
{
   power *= r;
   toSI = this->raise(toSI,r);
   this->symbol_to_the_power(r,1);
   if(symbol == "failed")
   {
     symbol = this->get_SI_symb();
   }

   return *this;
}

template <typename T>
Units<T> & Units<T>::operator/=(int r)
{
  power /= r;//check consistency of root
  this->symbol_to_the_power(r,-1);
  if(symbol == "failed")
  {
     symbol = this->get_SI_symb();
  }
  toSI = this->raise(toSI,1./(T)r);

  return *this;
}


template <typename T>
Units<T> Units<T>::operator+(const Units<T> & rhs) const
{
  return (Units<T>(*this) += rhs);
}

template <typename T>
Units<T> Units<T>::operator-(const Units<T> & rhs) const
{
  return (Units<T>(*this) -= rhs);
}

template <typename T>
Units<T> Units<T>::operator*(int r) const
{
  return (Units<T>(*this) *= r);
}

template <typename T>
Units<T> Units<T>::operator/(int r) const
{
  return (Units<T>(*this) /= r);
}
} //Antioch namespace

#endif
