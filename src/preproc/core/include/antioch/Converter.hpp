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
#ifndef _UNITS_CONVERTER_
#define _UNITS_CONVERTER_

//Antioch

//C++
#include <iostream>

namespace Antioch{

/*! \class Converter
 * \brief Class to deal with the conversion between units
 *
 * This class keeps track of the conversion factor between
 * the unit and the base unit. To be general, this is a couple
 * of double, one for the multiplication part and another one
 * for the translationnal part (mainly for temperature, but
 * it should not be used).
 */
class Converter{
      public:
/*! \brief Default constructor, defaults values are (1.,0.)*/
        Converter(double coef=1.,double trans=0.):a(coef),b(trans){}
/*! \brief Copy constructor, uses Converter & operator=(const Converter&)*/
        Converter(const Converter &rhs){*this = rhs;}
/*! \brief Default constructor*/
        ~Converter(){}

/*! \brief Multiplicative factor getter*/
        double geta() const {return a;}
/*! \brief Translationnal factor getter*/
        double getb() const {return b;}
/*! \brief Set the values to (1.,0.)*/
        void clear()        {a = 1.; b = 0.;}

/*! \brief << operator, to show the coef as a complex number (a,b)*/
        friend std::ostream &operator<< (std::ostream &out, const Converter & rhs)
        {
          out << "(" << rhs.geta() << "," << rhs.getb() << ")";
          return out;
        }

/*! \brief Comparison operator, equal if the two values are equal*/
        bool const operator== (const Converter &rhs) const;
/*! \brief Comparison operator, not equal is not "equal"*/
        bool const operator!=(const Converter &rhs) const {return !(*this == rhs);}
/*! \brief Assignement operator, equalize the two values*/
        Converter & operator=  (const Converter &rhs);
/*! \brief Multiplying operator
 *
 * The operations performed for the mutiplication of (a,b) and (c,d) are:
 * \f[
 *   a = a \times c
 * \f]
 * \f[
 *   b = \frac{b + d}{c}
 * \f]
 */
        Converter & operator*= (const Converter &rhs);
/*! \brief Dividing operator
 *
 * The operations performed for the division of (a,b) and (c,d) are:
 * \f[
 *   a = \frac{a}{c}
 * \f]
 * \f[
 *   b = \frac{b - d}{c}
 * \f]
 */
        Converter & operator/= (const Converter &rhs);
/*! \brief Adding operator
 *
 * The operations performed for the addition of (a,b) and (c,d) are:
 * \f[
 *   a = a \times c
 * \f]
 * \f[
 *   b = b + d
 * \f]
 */
        Converter & operator+= (const Converter &rhs);
/*! \brief Dividing operator, same operations as Converter &operator/=(const Converter&)*/
        Converter operator/    (const Converter &rhs) const;
/*! \brief Multiplying operator, same operations as Converter &operator*=(const Converter&)*/
        Converter operator*    (const Converter &rhs) const;
/*! \brief Multiplying operator with a double, only multiply the multiplicative value*/
        Converter & operator*= (double coef);
/*! \brief Dividing operator with a double, only divide the multiplicative value*/
        Converter & operator/= (double coef);
/*! \brief Multiplying operator with a double, only multiply the multiplicative value*/
        Converter operator*    (double coef)   const;
/*! \brief Dividing operator with a double, only divide the multiplicative value*/
        Converter operator/    (double coef)   const;

      private:
/*! \brief Two double for a multiplicative and translationnal part of a coefficient*/
        double a,b;

};
}

#endif
