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
#ifndef ANTIOCH_CONVERTER_H
#define ANTIOCH_CONVERTER_H

//Antioch
#include "antioch/metaprogramming.h"

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
template <typename T = double>
class Converter{
      public:
/*! \brief Default constructor, defaults values are (1.,0.)*/
        Converter();
/*! \brief Constructor*/
        Converter(const T &coef, const T &trans);
/*! \brief Copy constructor, uses Converter & operator=(const Converter&)*/
        template <typename P>
        Converter(const Converter<P> &rhs){*this = rhs;}
/*! \brief Default constructor*/
        ~Converter(){}

/*! \brief Multiplicative factor getter*/
        template <typename P = T>
        const P geta() const {return (P)a;} //only scalar type
/*! \brief Translationnal factor getter*/
        template <typename P = T>
        const P getb() const {return (P)b;} //only scalar type
/*! \brief Set the values to (1.,0.)*/
        void clear()        {a = 1.L; Antioch::set_zero(b);}

/*! \brief << operator, to show the coef as a complex number (a,b)*/
        friend std::ostream &operator<< (std::ostream &out, const Converter<T> & rhs)
        {
          out << "(" << rhs.geta() << "," << rhs.getb() << ")";
          return out;
        }

/*! \brief Comparison operator, equal if the two values are equal*/
        template <typename P>
        bool operator== (const Converter<P> &rhs) const;
/*! \brief Comparison operator, not equal is not "equal"*/
        template <typename P>
        bool operator!=(const Converter<P> &rhs) const {return !(*this == rhs);}
/*! \brief Assignement operator, equalize the two values*/
        template <typename P>
        Converter & operator=  (const Converter<P> &rhs);
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
        template <typename P>
        Converter & operator*= (const Converter<P> &rhs);
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
        template <typename P>
        Converter & operator/= (const Converter<P> &rhs);
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
        template <typename P>
        Converter & operator+= (const Converter<P> &rhs);
/*! \brief Dividing operator, same operations as Converter &operator/=(const Converter&)*/
        template <typename P>
        Converter operator/    (const Converter<P> &rhs) const;
/*! \brief Multiplying operator, same operations as Converter &operator*=(const Converter&)*/
        template <typename P>
        Converter operator*    (const Converter<P> &rhs) const;
/*! \brief Multiplying operator with a scalar, only multiply the multiplicative value*/
        template <typename P>
        Converter & operator*= (const P &coef);
/*! \brief Dividing operator with a double, only divide the multiplicative value*/
        template <typename P>
        Converter & operator/= (const P &coef);
/*! \brief Multiplying operator with a scalar, only multiply the multiplicative value*/
        template <typename P>
        Converter operator*    (const P &coef)   const;
/*! \brief Dividing operator with a scalar, only divide the multiplicative value*/
        template <typename P>
        Converter operator/    (const P &coef)   const;

      private:
/*! \brief Two double for a multiplicative and translationnal part of a coefficient*/
        T a;
        T b;

};

template<typename T>
inline
Converter<T>::Converter()
{
  clear();
  return;
}

template<typename T>
inline
Converter<T>::Converter(const T &coef, const T &tran)
{
  a = Antioch::constant_clone(a,coef);
  b = Antioch::constant_clone(b,tran);
  return;
}


template<typename T>
template<typename P>
inline
Converter<T> & Converter<T>::operator= (const Converter<P> &rhs)
{
 // if(this == static_cast<Converter<T>*>(&rhs)){return *this;}
  a = Antioch::constant_clone(a,rhs.geta());
  b = Antioch::constant_clone(b,rhs.getb());
  return *this;
}

template<typename T>
template<typename P>
inline
Converter<T> & Converter<T>::operator*= (const P &coef)
{
  a *= coef;
  return *this;
}

template<typename T>
template<typename P>
inline
Converter<T> & Converter<T>::operator/= (const P &coef)
{
  a /= coef;
  return *this;
}

template<typename T>
template<typename P>
inline
Converter<T> & Converter<T>::operator*= (const Converter<P> &rhs)
{
  a *= rhs.geta();
  b  = (b + rhs.getb())/rhs.geta();
  return *this;
}

template<typename T>
template<typename P>
inline
Converter<T> & Converter<T>::operator/= (const Converter<P> &rhs)
{
  a /= rhs.geta();
  b  = (b - rhs.getb())/rhs.geta();
  return *this;
}

template<typename T>
template<typename P>
inline
Converter<T> & Converter<T>::operator+= (const Converter<P> &rhs)
{
  a *= rhs.geta();
  b += rhs.getb();
  return *this;
}

template<typename T>
template<typename P>
inline
Converter<T> Converter<T>::operator*  (const P &coef) const
{
  return Converter<T>(a*coef,b);
}

template<typename T>
template<typename P>
inline
Converter<T> Converter<T>::operator/  (const P &coef) const
{
  return Converter<T>(a/coef,b);
}

template<typename T>
template<typename P>
inline
Converter<T> Converter<T>::operator/  (const Converter<P> &rhs) const
{
  return Converter(a/rhs.geta(),(b-rhs.getb())/rhs.geta());
}

template<typename T>
template<typename P>
inline
Converter<T> Converter<T>::operator*  (const Converter<P> &rhs) const
{
  return Converter<T>(a * rhs.geta() , (b + rhs.getb())/rhs.geta());
}
}// end namespace

#endif
