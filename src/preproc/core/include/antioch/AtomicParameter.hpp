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

#ifndef _ATOMIC_PARAMETER_
#define _ATOMIC_PARAMETER_

//Antioch

//C++
#include <string>
#include <vector>

namespace Antioch{

/*! \class AtomicParameter
 *  \brief An atomic parameter is the least you can have as a physical parameter
 *
 * An AtomicParameter stores:
 *   - a name
 *   - associated value(s)
 *
 * It does not do anything, only stores
 */
class AtomicParameter{
  public:
/*! \brief Copy constructor, use the AtomicParameter &operator=(const AtomicParameter&)*/
   AtomicParameter(const AtomicParameter &ToBeCopied){*this = ToBeCopied;}
/*! \brief Default constructor, a std::string is given for the name of the parameter*/
   AtomicParameter(std::string namePar = std::string()):name(namePar){}
/*! \brief Building constructor, a std::string and a std::vector of doubles for the values*/
   AtomicParameter(std::string namePar,std::vector<double> valuesPar):name(namePar),values(valuesPar){}
/*! \brief Building constructor, a std::string and a double for a value*/
   AtomicParameter(std::string namePar,double firstValuePar):name(namePar){values.push_back(firstValuePar);}
/*! \brief Building constructor,  only a double for a value*/
   AtomicParameter(double firstValuePar):name("scalar"){values.push_back(firstValuePar);}
/*! \brief Default destructor*/
   ~AtomicParameter(){}

/*! \brief String name getter*/
   const std::string getName()             const {return name;}
/*! \brief String name setter*/
   void setName(const std::string &namePar)      {name = namePar;}
/*! \brief Name clearer*/
   void clearName()                              {name.clear();}

/*! \brief Value getter*/
   double getValue(int i)                  const;
/*! \brief Overload of double getValue(int)*/
   double getValue(unsigned int i)         const {return getValue((int)i);}
/*!\brief First value getter*/
   double front()                          const {return values.front();}
/*!\brief Last value getter*/
   double back()                           const {return values.back();}
/*! \brief Value setter
 *
 * By default, if not provided, the value 0 is set (or added if
 * need be), so it eases constants handling.
 * */
   void setValue(double val,int ival = 0);
/*! \brief Overload of void setValue(double, int)*/
   void setValue(double val,unsigned int ival)         {setValue(val,(int)ival);}
/*! \brief Value adder*/
   void addValue(double val)                           {values.push_back(val);}

/*! \brief Values setter*/
   void setValues(const std::vector<double> &vals)     {values = vals;}
/*! \brief Values adder*/
   void addValues(const std::vector<double> &vals);
/*! \brief Values getter*/
   const std::vector<double> getValues()         const {return values;}

/*\brief Duplicate a value a given number of times*/
   void duplicateVal(int ival, int ntimes);
/*\brief Duplicate a value a given number of times, unsigned int overload*/
   void duplicateVal(unsigned int ival, unsigned int ntimes) {duplicateVal((int)ival,(int)ntimes);}

/*! \brief Number of values getter*/
   const int sizeVal()                      const {return values.size();}
/*! \brief Values clearer*/
   void clearValues()                             {values.clear();}
/*! \brief Values resizer*/
   void resizeValues(unsigned int newSize, double val) {values.resize(newSize,val);}
/*! \brief Overload of void resizeValues(unsigned int)*/
   void resizeValues(int newSize, double val)          {resizeValues((unsigned int)newSize,val);}
/*! \brief Values reserver*/
   void reserveValues(unsigned int newSize)       {values.reserve(newSize);}
/*! \brief Overload of void reserveValues(unsigned int)*/
   void reserveValues(int newSize)               {reserveValues((unsigned int)newSize);}

/*! \brief Boolean equal operator
 *
 * Two AtomicParameter are equal if and only if
 * their names are the same and the values are
 * the same, in the same order.
 * */
   const bool operator==(const AtomicParameter &rhs) const;
/*! \brief Not const bool operator==(const AtomicParameter&)*/
   const bool operator!=(const AtomicParameter &rhs) const;
/*! \brief Assignement operator, the name and the values are copied*/
   AtomicParameter & operator =(const AtomicParameter &rhs);
/*! \brief Multiplying operator
 *
 * Some checks are necessary before multiplying the values:
 *   - there should be values
 *   - the parameters have the same number of values, or
 *   one parameter at least should have only one value
 */
   AtomicParameter & operator*=(const AtomicParameter &rhs);
/*! \brief Multiplying operator, uses AtomicParameter &operator*=(const AtomicParameter&)*/
   AtomicParameter   operator* (const AtomicParameter &rhs) const;
/*! \brief Dividing operator, dividing equivalent to AtomicParameter &operator*=(const AtomicParameter&)*/
   AtomicParameter & operator/=(const AtomicParameter &rhs);
/*! \brief Dividing operator, uses AtomicParameter &operator/=(const AtomicParameter&)*/
   AtomicParameter   operator/ (const AtomicParameter &rhs) const;
/*! \brief Adding operator, adding equivalent to AtomicParameter &operator*=(const AtomicParameter&)*/
   AtomicParameter & operator+=(const AtomicParameter &rhs);
/*! \brief Adding operator, uses AtomicParameter &operator+=(const AtomicParameter&)*/
   AtomicParameter   operator+ (const AtomicParameter &rhs) const;
/*! \brief Substracting operator, substracting equivalent to AtomicParameter &operator*=(const AtomicParameter&)*/
   AtomicParameter & operator-=(const AtomicParameter &rhs);
/*! \brief Substracting operator, uses AtomicParameter &operator-=(const AtomicParameter&)*/
   AtomicParameter   operator- (const AtomicParameter &rhs) const;
/*! \brief Multiplying operator with a double, simple multiplication*/
   AtomicParameter & operator*=(double r);
/*! \brief Multiplying operator with a double, 
 * uses AtomicParameter &operator*=(double)*/
   AtomicParameter   operator* (double r) const;
/*! \brief Dividing operator with a double, this is a multiplication by the inverse of the double*/
   AtomicParameter & operator/=(double r);
/*! \brief Dividing operator with a double, 
 * uses AtomicParameter &operator/=(double)*/
   AtomicParameter   operator/ (double r) const;
/*! \brief Adding operator with a double, simple addition*/
   AtomicParameter & operator+=(double r);
/*! \brief Adding operator with a double, 
 * uses AtomicParameter &operator+=(double)*/
   AtomicParameter   operator+ (double r) const;
/*! \brief Substracting operator with a double, simple substraction*/
   AtomicParameter & operator-=(double r);
/*! \brief Substracting operator with a double, 
 * uses AtomicParameter &operator-=(double)*/
   AtomicParameter   operator- (double r) const;
/*! \brief Inversing the values */
   AtomicParameter   operator- () const;

/*! \brief Alternative call to operator*/
   void equalize(const AtomicParameter &rhs)  {*this  = rhs;}
/*! \brief Alternative call to operator*/
   void add(const AtomicParameter &rhs)       {*this += rhs;}
/*! \brief Alternative call to operator*/
   void substract(const AtomicParameter &rhs) {*this -= rhs;}
/*! \brief Alternative call to operator*/
   void multiply(const AtomicParameter &rhs)  {*this *= rhs;}
/*! \brief Alternative call to operator*/
   void divide(const AtomicParameter &rhs)    {*this /= rhs;}
/*! \brief Alternative call to operator*/
   void add(double r)                         {*this += r;}
/*! \brief Alternative call to operator*/
   void substract(double r)                   {*this -= r;}
/*! \brief Alternative call to operator*/
   void multiply(double r)                    {*this *= r;}
/*! \brief Alternative call to operator*/
   void divide(double r)                      {*this /= r;}

/*! \brief Clear functions, clear the values and name*/
   void clear()      {clearValues();clearName();}

  private:

/*! \brief Intermediary method 
 *
 * It is used to check the parameters for
 * computation. Checks that:
 *   - one parameter is a constant, thus only one value, or,
 *   - the two parameters have the same number of values
 *
 * It returns the maximum numbers of values for computation
 */
   int testArithmOperation(const std::string &method,const AtomicParameter &rhs, const std::string &op, const std::string &sym) const;

/*! \brief Intermediary method 
 *
 * It is useful for error handling, produce the long std::string
 * in case of error in an operator
 */
   const std::string computationError(std::string ope, std::string op) const;

/*! \brief A std::string for a name*/
   std::string name;
/*! \brief A std::vector<double> for the values*/
   std::vector<double> values;
};

}

#endif
