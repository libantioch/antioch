//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
#ifndef _UNITS_MANAGER_
#define _UNITS_MANAGER_

//Antioch
#include "antioch/Error.hpp"
#include "antioch/SIPrefix.hpp"
#include "antioch/InSI.hpp"
#include "antioch/Converter.hpp"

//C++
#include <sstream>
#include <iomanip>

namespace Antioch{

/*!\file Units.hpp
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
 * At each unit parsing level, the unit is parsed between a prefixe, a symbol and a power.
 * The symbol "cm-3" will be parsed as a prefixe "c", a symbol "m" and a power "-3". To this
 * parsing, the symbol before the unit is considered. If there is none or it is a dot, nothing
 * to do, but if this is a slash, the power has to be multiplied by -1. Thus
 * the Converter object will have the corresponding SIPrefixes values (Prefixe[6] is this case,
 * which is ("c",1e-2)), and the InSI power vector will add the value "-3" to the corresponding
 * power (0 in this case) if the character before the unit is a dot or does not exist,
 * the value 3 if this character is a slash. The correspondancies between the index in the power vector and
 * the unit is done explicitely.
 *
 * A special treatment is reserved for the kg unit, as it contains a prefixe.
 *
 * This method supports a wide variety of symbol writing, like the use of parenthesises. The
 * used symbol is a fully developped one. 
 */
class Units{
      public:
/*!\brief Copy constructor, uses Units &operator=(const Units&)*/
        Units(const Units &rhs){*this = rhs;developSymbol(symbol);}
/*!\brief Fully descriptive constructor, mainly for internal purposes*/
        Units(std::string sym,std::string na,
                double conva,double convb,
                int mi,int kgi=0, int si=0, int Ai=0, int Ki=0, int moli=0, int cdi=0, int radi=0);//constructors for known Units
/*!\brief User-friendly constructor, from a symbol and an optionnal name*/
        Units(std::string sym,std::string na="");  //constructors when automatic conversion
/*!\brief Constructor for a user imposed conversion*/
        Units(std::string sym,Converter conv, std::string na="");//constructor for given conversion
/*!\brief Default constructor*/
        Units(){}

/*!\brief Default destructor*/
        ~Units(){}

/*! \brief Comparison operator
 *
 * Compares the power vector and the coefficient, regardless
 * of the symbol and the name.
 */
        bool const operator== (const Units &rhs) const;
/*! \brief Comparison operator, not equal is not "equal"*/
        bool const operator!= (const Units &rhs) const {return !(*this == rhs);}
/*!\brief Assignement operator
 *
 * It assigns the symbol std::string,
 * the power vector and coefficient are filled.
 */
        Units &operator= (const Units & rhs);
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
        Units &operator+=(const Units & rhs);
/*!\brief Substracting operator
 *
 * Substracting units occurs when dividing physical quantities.
 *
 * Substracting is close to adding (see Units &operator+=(const Units&)),
 * with little differences. The adding of the name is done
 * with the sequence blank-slash-blank (" / ") between the names.
 * The symbol needs a bit more processing: if the added symbol is
 * not empty, is it put between parenthesises and added with a slash
 * in between. Then the void developSymbol() method is called to
 * suppress the parenthesises and obtain a correct symbol.
 *
 * The coefficient is divided and
 * the power vector substracted.
 */
        Units &operator-=(const Units & rhs);
/*!\brief Adding operator, uses Units &operator+=(const Units&)*/
        Units  operator+ (const Units & rhs) const;
/*!\brief Substracting operator, uses Units &operator-=(const Units&)*/
        Units  operator- (const Units & rhs) const;
/*!\brief Multiplying operator by an integer
 *
 * This operator is needed for power elevation operation of a physical quantity.
 *
 * The operations needed are as follow:
 *   - parsing the symbol so all the power are multiplied by the integer,
 *   done by the method void symbolToThePower(int,std::string);
 *   - multiplying the power vector by the integer;
 *   - elevating the coefficient to the power, done by the method
 *   Converter raise(const Converter&, int).
 *
 * Elevating to the power with a double is not supported, and the physical sense
 * of it anyway is not obvious...
 */
        Units &operator*=(int r);
/*!\brief Multiplying operator by an integer, uses Units &operator*=(int)*/
        Units  operator* (int r)             const;
/*!\brief dividing operator by an integer
 *
 * This operator is needed for root operation of a physical quantity.
 *
 * The operations needed are as follow:
 *   - dividing the power vector by the integer, test that the resulting
 *     power is an integer;
 *   - parsing the symbol so all the power are divided by the integer,
 *     done by the method void symbolToThePower(int,std::string),
 *     if the symbol contains unit symbol that gives non integer
 *     resulting power, then the SI symbol is called;
 *   - dividing the coefficient to the power, done by the method
 *     Converter raise(const Converter&, double).
 */
        Units &operator/=(int r);
/*!\brief Dividing operator by an integer, uses Units &operator/=(int)*/
        Units  operator/ (int r)             const;

/*!\brief Alternative call to Units& operator=(const Units &) */
        void equalize(const Units &rhs)      {*this  = rhs;}
/*!\brief Alternative call to Units& operator=(const Units &) */
        void equalize(Units *rhs)            {*this  = *rhs;}
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
        void equalize(std::string unit)       {*this  = Units(unit);}
/*!\brief Alternative call to Units& operator+=(const Units &) */
        void add(const Units &rhs)            {*this += rhs;}
/*!\brief Alternative call to Units& operator+=(const Units &) */
        void add(Units *rhs)                  {*this += *rhs;}
/*!\brief Alternative call to Units& operator+=(const Units &) 
 *
 * Same remarks than void equalize(std::string).
 */
        void add(std::string unit)            {*this += Units(unit);}
/*!\brief Alternative call to Units& operator-=(const Units &)*/
        void substract(const Units &rhs)      {*this -= rhs;}
/*!\brief Alternative call to Units& operator-=(const Units &)*/
        void substract(Units *rhs)            {*this -= *rhs;}
/*!\brief Alternative call to Units& operator-=(const Units &) 
 *
 * Same remarks than void equalize(std::string).
 */
        void substract(std::string unit)      {*this -= Units(unit);}
/*!\brief Alternative method to Units& operator/=(int)*/
        void root(int r)                      {*this /= r;}

/*! \brief Units setter, sets the unit to the given symbol and name*/
        void setUnit(const std::string &sym, std::string na);
/*! \brief Units setter, sets the unit to the given symbol and name*/
        void setUnit(const std::string &sym)       {setUnit(sym,"");}
/*! \brief Name setter*/
        void setName(const std::string &na)        {name = na;}
/*! \brief Symbol setter*/
        void setSymbol(const std::string &symb)    {symbol = symb;}
/*! \brief Coefficient setter*/
        void setSIConverter(const Converter &conv) {toSI = conv;}
/*! \brief Power std::vector setter*/
        void setPowerToSI(const InSI &mat)         {power = mat;}

/*!\brief Homogenity testing with another Units.
 *
 * Testing the homogeneity between units is different from testing
 * equality. J and cal are homogeneous but not equal for instance.
 * Thus this method merely compares the power vectors, if they're equal
 * the units are homogeneous, if not, they're not.
 */
        const bool isHomogeneous(const Units &rhs)    const {return (power == rhs.getPower());}
/*!\brief Homogenity testing with a std::string.
 *
 * Testing homogeneity with a std::string needs to build
 * a Units object and testing homogeneity with it. In the
 * case of an empty std::string, instead of building a Units object
 * the boolean isUnited() is tested.
 */
        const bool isHomogeneous(std::string target) const;
/*!\brief Test if the unit is non empty*/
        const bool isUnited() const {return !power.empty();}


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
 * given above. FactorToSomeUnit(unit\f$_2\f$) will provide the factor \f$\frac{f_1}{f_2}\f$
 * while TranslatorToSomeUnit(unit\f$_2\f$) will provide the translationnal term \f$\frac{t_1 - t_2}{f_2}\f$.
 */
        double FactorToSomeUnit(const Units & target)  const;
/*!\brief Calculates the factor to any given unit, uses FactorToSomeUnit(const Units&).*/
        double FactorToSomeUnit(std::string target)         const {return FactorToSomeUnit(Units(target));}
/*!\brief Calculates the translator to any given unit, see FactorToSomeUnit(const Units&) for the equations.*/
        double TranslatorToSomeUnit(const Units & target)  const;
/*!\brief Calculates the translator to any given unit, uses TranslatorToSomeUnit(const Units&).*/
        double TranslatorToSomeUnit(std::string target)  const {return TranslatorToSomeUnit(Units(target));}

/*! \brief String name getter*/
        const std::string getName()       const {return name;}
/*! \brief String symbol getter*/
        const std::string getSymbol()     const {return symbol;}
/*! \brief Coefficient getter*/
        const Converter getSIcoef()       const {return toSI;}
/*! \brief Multiplicative coefficient getter*/
        double getSIFactor()              const {return toSI.geta();}
/*! \brief Translationnal coefficient getter*/
        double getSITranslator()          const {return toSI.getb();}
/*! \brief Power std::vector getter*/
        const InSI getPower()             const {return power;}
/*! \brief Corresponding SI symbol getter
 *
 * This method returns the symbol that correspond to the power vector.
 * For instance a unit J/s will be returned as m2.kg.s-3
 */
        const std::string getSISymb()     const;
/*! \brief Power of the asked SI unit
 *
 * This method returns the power to the 
 * asked unit SI.
 */
        int getSIPower(const std::string &SIask) const;
/*! \brief Corresponding SI symbol getter
 *
 * This method returns the symbol that of the SI unit (so a
 * coefficient of (1.,0.)). If none is found, it returns
 * const std::string getSISymb().
 * For instance a unit J/s will be returned as W.
 *
 * This method is not able to combine different SI units.
 * The unit J.s will not be recognized as J.s, only as
 * m2.kg.s-1
 */
        const std::string getSIConvenientSymb() const;
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
        const std::string harmonizedSymbol(const std::string &input = "")    const {return manipulateSymbol(input,false);}
/*! \brief Contract the symbol when possible.
 *
 * The contraction will aggregate only exactly equal
 * units.
 */
        const std::string contractedSymbol(const std::string &input = "")    const {return manipulateSymbol(input,true);}

/*! \brief Clear the unit.*/
        void clear(){symbol.clear();name.clear();toSI.clear();power.clear();}

/*! \brief showAll().*/
        void showAll(std::ostream &out = std::cout);


      private:

/*! \brief Root method for contracting or harmonizing the symbol
 *
 *  The symbol MUST be developped prior to be used in this method.
 */
        const std::string manipulateSymbol(std::string input, bool contract) const;


/*!\brief This method fills the power vector and will calculate the coefficient
 * if the bool doConv is set to true.
 *
 * The method scans the symbol, and calls bool parseSingleUnit(int,std::string,bool)
 * for each unit found. The character in front of the unit ('.' or '/') is taken
 * into account into the first integer (resp. 1 or -1).
 * This is this latter method that actually calculates the coefficients and
 * power vector.
 */
        void fillInPower(bool doConv);
/*!\brief Calculates the corresponding coefficient and power of the given unit.
 *
 * \param int signe, 1 or -1. Image of the sign before the unit '.' or '/'.
 * \param std::string unit. The unit to be parsed.
 * \param bool doConv. Switch to perform or not the coefficient calculation.
 *
 * The method works in two times:
 *   - first parsing the symbol into prefixe + symbol + power
 *   - second update the power vector and, if doConv is true, the coefficient.
 * 
 * The parsing uses the methods int parsePower(std::string,int&) and parsePrefixeUnit(int&,int&,std::string).
 * The idea is to scan the knownUnits[] and Prefixes[] array and find the corresponding unit
 * and prefixe. The method int parsePower returns the power, while
 * the integer iPre and iUnit characterize the prefixe and unit, with a value of -1 if
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
 * coefficient of the prefixe, raise to the power by the parsed power, updated by the signe:
 * \f[
 *      (c) = (c) \times \left( pr_{iPre} \times (c)_{iUnit}\right)^{s_u \times p_u}
 * \f]
 * with \f$(c)\f$ being the updated complex coefficient, \f$(c)_{iUnit}\f$ the complex
 * coefficient of the found unit, \f$pr_{iPre}\f$ the coefficient of the
 * found Prefixes (if found, if not, it is unity), 
 * \f$s_u\f$ the sign of the parsed unit and \f$p_u\f$ the
 * power of the parsed unit.
 */
        bool parseSingleUnit(int signe, std::string unit, bool doConv);
/*!\brief Unit parser, parse a std::string of type Prefixes[] knownUnits[].
 *
 * \param int &iUnit. Output, indice of the found unit, if not found, set to -1.
 * \param int &iPre. Output, indice of the found prefixes, if not found, set to -1.
 * \param std::string unit. Input, std::string to be parsed.
 *
 * This method will parse a std::string, assuming the std::string is of the form
 * Prefixes[]knownUnits[]. No power should be present (see int parsePower(std::string,int&)
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
 * (including no prefixe) are different.
 */
        void parsePrefixeUnit(int &iUnit,int &iPre,std::string unit)  const;
/*!\brief Scanner of knownUnits[], sends back the corresponding iUnit, -1 if not found.*/
        int indexUnit(std::string unit)                               const;
/*!brief parser for std::string of type "some_unit power_associated".
 *
 * \param std::string unit. Input, std::string to be parsed.
 * \param int &nc. Output, number of charaters of the power.
 * \return the power associated to the unit.
 */
        int parsePower(std::string unit, int &nc)                    const;
/*!\brief Small method to check if a character is an alpha-numerical.
 *
 * This method uses the values of the ascii table, which basically
 * means that is checks if the ascii integer associated to the
 * character is in the range [48,57].
 */
        bool isNumber(char c)                                   const {return ((int)c <= 57 && (int)c >= 48);} //ascii table

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
        Converter raise(const Converter &tbm,int power)       const;
/*!\brief Overload to simplify division
 *
 * Another method, to secure multiplication. Division is thus
 * put apart, to secure against non integer power.
 */
        Converter raise(const Converter &tbm,double power)       const;
/*!\brief Supress the parenthesises
 *
 * Two types of parenthesises are to be treated:
 *   - serial parenthesises ().()
 *   - intricated parenthesises (())
 *
 * The idea is to treat every serial pairs one by one.
 * Thus every pair found is
 * sent back to developSymbol(std::string&). Once the
 * we have treated the inside string by this recursive method,
 * we treat the string:
 *   - if the character before the opening parenthesis is '.',
 *      we just supress the parenthesises;
 *   - if it is '/', we reversePowerSymbol(std::string &) the string
 */
        void developSymbol(std::string &subsymbol);
/*!\brief Treat part of symbol depending on parenthesises
 *
 * Inverse the power sign determinant ('.' or '/'). Note that
 * the first power will not be affected. Ex:
 * J.mol-1/K will be changed to J/mol-1.K, leaving
 * J unchanged.
 */
        void reversePowerSymbol(std::string &subsymbol);
/*!\brief Small method that checks if the symbol is within a unit std::string.
 *
 * This method checks if the symbol is not a number (bool isNumber(char)),
 * and not any of the following character: '/', '.' and '-'. Then this
 * symbol is consider to belong to a unit symbol.
 */
        const bool isInSymb(char c)               const;
/*!\brief Number of dimensions of the unit in the SI basis.*/
        const int nDimensionOfUnits()             const;
/*!\brief Small method to add the SI symbol and power to a std::string.
 *
 * This method is useful in the std::string const Units::getSISymb() const
 * method.
 */
        const std::string addSI(int pow,std::string SIsymb) const;
/*!\brief Rewrite the symbol to put it at a given power.
 *
 * This method scans the symbol, parse the units, and changes
 * the power by the power multiplied or divided by the given integer.
 * It uses the method getIntegerPower(int,int,std::string) to ensure
 * integer power.
 */
        void symbolToThePower(int r, const std::string &key);
/*!\brief Multiply or divide a power
 *
 * It ensures that the result of the operation performed
 * is an integer. Natural in case of multiplication, to be
 * checked in case of division. It returns 0 if it is
 * not correct, this \e has to be checked by the wrapping
 * method.
 * The key is to choose between division or multiplication.
 * Key is:
 *   - ``multiplication''
 *   - ``division''
 */
        int getIntegerPower(int unit,int r, const std::string &key);
/*!\brief Fill the name if empty and name obvious
 *
 * Usually we give only the symbol, even for simple unit.
 * Thus this method will fill the name if the name is obvious,
 * i.e. if the dimension of the power vector is 1. The
 * suffixe is given by the factor value.
 */
        void checkIfName(std::string &nameOut, const Units &un) const;

/*!\brief Strings for the symbol.*/
        std::string symbol;
/*!\brief Strings for the name.*/
        std::string name;
/*!\brief A Converter for the coefficient.*/
        Converter toSI;
/*!\brief An InSI for the power vector.*/
        InSI power;
};


/*!\brief Check if there is no redundancy in the known unit.
 *
 * It tests if each couple prefixe-unit is unique. This method
 * should be run each time there's a modification of KnownUnits[]
 * and/or Prefixes[].
 */
void AreDefsOk();

/*!\brief List all known Units, public function*/
std::string allKnownUnits();

/*!\brief List all known Units that are homogeneous to the given one, public function*/
std::string allKnownHomogeneousUnits(std::string target);

/*!\brief List all known prefixes, public function*/
std::string allKnownPrefixes();
}
#endif
