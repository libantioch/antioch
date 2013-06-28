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
#ifndef _PHYMATHS_FUNCTIONS_
#define _PHYMATHS_FUNCTIONS_

//Antioch
#include "antioch/ParameterPhy.hpp"

//C++

namespace Antioch{
/*!\file Phymath.hpp
 * \brief Classical maths functions with ParameterPhy objects
 * \author \me
 *
 * This file regroups public functions to be used
 * with ParameterPhy objects. This file regroups 
 * the following functions:
 *   - cos
 *   - sin
 *   - tan
 *   - acos
 *   - asin
 *   - atan
 *   - exp
 *   - log
 *   - log10
 *   - pow
 *   - fabs
 *   - erf
 *   - sqrt
 *
 * All those functions check for unit consistency
 * and perform uncertainty calculations, following the
 * GUM recommandations:
 * \f[
 *      u_c^2 = \sum_i \left(\frac{\partial f}{\partial x_i}\right)^2 u_i^2
 * \f]
 * with \f$u_i\f$ being the standard deviation of variable \f$x_i\f$
 * The uncertainty must be changed to CORE_UNCERTAINTY_TYPE_ABSOLUTE
 * if it is not the case.
 *
 * There is also constant parameters handling that must be performed.
 * A constant is a ParameterPhy with exactly one value. In the case
 * there is a constant to be combined with a non constant, proper
 * indices handling is performed.
 *
 * This is done with generic public method. Those methods ensure
 * standardized ParameterPhy calculations, for any function needed
 * not in this list. A function calculation is divided in three parts:
 *   - unit checks and coefficient definition
 *   - values calculations
 *   - uncertainty values calculations
 *
 */

/*!\brief Public method to check for unit consistency and coefficient definition
 *
 * The coefficient is necessary a multiplicative coefficient.
 * This coefficient is one, then the method
 * checks if the unit and unitRef are the same,
 * if no, it checks for homogeneity. If they are not
 * homogeneous, it sends back an error, if they are,
 * it returns the coefficient needed to pass from the unit
 * to unitRef.
 */
double checkFunctionUnit(const std::string unitRef, 
                         const std::string methodCalling, 
                         const std::string errmessage,
                         const ParameterPhy &rhs);


/*!\brief This is a wrapper function
 *
 * It performs \f$y_i = f(x_i \cdot c)\f$ for every
 * \f$y_i\f$ value of the given ParameterPhy \f$\mathbf{Y}\f$ with
 * the given function \f$f\f$ applied to every \f$x_i\f$ values
 * of given ParameterPhy \f$\mathbf{X}\f$, multiplied by the given coefficient
 * \f$c\f$.
 */
void calculateFunctionValues(ParameterPhy &y, 
                             const ParameterPhy &x, 
                             double coef, 
                             double function(double ));

/*!\brief This is a wrapper function
 *
 * It performs 
 *  \f$u_{y_i}^2 = \left(\frac{\mathrm{d}f}{\mathrm{d}x} (x_i \cdot c)\right)^2u_{x_i}^2\f$ 
 * for every \f$y_i\f$ value of the given ParameterPhy \f$\mathbf{Y}\f$ with
 * the given function \f$\frac{\mathrm{d}f}{\mathrm{d}x}\f$ applied to every \f$x_i\f$ values
 * of given ParameterPhy \f$\mathbf{X}\f$, multiplied by the given coefficient
 * \f$c\f$.
 * The uncertainty type is managed as follow:
 *   - CORE_UNCERTAINTY_TYPE_NONE and CORE_UNCERTAINTY_TYPE_UNKNOWN will 
 *     set the same uncertainty type to the ParameterPhy \f$\mathbf{Y}\f$ and
 *     do nothing
 *   - CORE_UNCERTAINTY_TYPE_RELATIVE will perform:
 *      \f[
 *           u_{y_i}^2 = \frac{\mathrm{d}f}{\mathrm{d}x}(x_i \cdot c) \cdot 
 *                        \frac{\mathrm{d}f}{\mathrm{d}x}(x_i \cdot c) \cdot 
 *                        u_{x_i} \cdot c \cdot x_i \cdot
 *                        u_{x_i} \cdot c \cdot x_i 
 *      \f]
 *   - CORE_UNCERTAINTY_TYPE_ABSOLUTE  will perform:
 *      \f[
 *           u_{y_i}^2 = \frac{\mathrm{d}f}{\mathrm{d}x}(x_i \cdot c) \cdot 
 *                        \frac{\mathrm{d}f}{\mathrm{d}x}(x_i \cdot c) \cdot 
 *                        u_{x_i}^2 \cdot c \cdot c
 *      \f]
 *      
 *  Remember that this is the variance that is stored in case of CORE_UNCERTAINTY_TYPE_ABSOLUTE 
 *  and the relative value in the CORE_UNCERTAINTY_TYPE_RELATIVE case. The uncertainty type
 *  of \f$\mathbf{Y}\f$ is always CORE_UNCERTAINTY_TYPE_ABSOLUTE.
 *
 *  It is noticeable that giving \f$f(x)\f$ or \f$-f(x)\f$ is the same.
 */
void calculateFunctionDValues(ParameterPhy &y, 
                             const ParameterPhy &x, 
                             double coef, 
                             double function(double ));


/* \brief Derived function for tangente
 *
 * This is \f$\frac{\mathrm{d}tan}{\mathrm{d}x} = \frac{1}{\cos^2(x)}\f$
 */
double dtandx(double x);

/* \brief Derived function for arccosinus
 *
 * This is \f$\frac{\mathrm{d}acos}{\mathrm{d}x} = -\frac{1}{\sqrt{1 - x^2}}\f$
 */
double dacosdx(double x);

/* \brief Derived function for arcsinus
 *
 * This is \f$\frac{\mathrm{d}acos}{\mathrm{d}x} = \frac{1}{\sqrt{1 - x^2}}\f$
 */
double dasindx(double x);

/* \brief Derived function for natural logarithm
 *
 * This is \f$\frac{\mathrm{d}\ln}{\mathrm{d}x} = \frac{1}{x}}\f$
 */
double dlogdx(double x);

/* \brief Derived function for base 10 logarithm
 *
 * This is \f$\frac{\mathrm{d}\ln}{\mathrm{d}x} = \frac{1}{x \ln(10)}}\f$
 */
double dlog10dx(double x);

/*!\brief Error function
 *
 * The approximation from Abramowitz and Stegun, Eq. 7.1.26
 * \f[
 *    0 \leq z < \infty \qquad  \mathrm{erf}(z) = 1 - (a_1 t + a_2 r^2 + a_3 t^3 + a_4 t^4 + a_5 t^5)e^{-z^2} + \epsilon(z)
 * \f]
 * with
 *   - \f$t = \frac{1}{1 + p z}\f$
 *   - \f$p = 0.3275911\f$
 *   - \f$a_1 = 0.254829592\f$
 *   - \f$a_2 = -0.284496736\f$
 *   - \f$a_3 = 1.421413741\f$
 *   - \f$a_4 = -1.453152027\f$
 *   - \f$a_5 = 1.061405429\f$
 *   - \f$|\epsilon(z)| \leq 1.5\,10^{-7}\f$
 *
 * For negative values, we used the fact that 
 * \f$\mathrm{erf}(-z) = -\mathrm{erf}(z)\f$
 */
double erf(double z);

/*! \brief Derived function for error function
 *
 * This is \f$\frac{\mathrm{d}\mathrm{erf}}{\mathrm{d}x} = \frac{2}{\sqrt{\pi}} \exp(-x^2)\f$
 */
double derfdx(double x);

/*! \brief Derived of root square function
 *
 * This is \f$\frac{1}{2\sqrt{x}}\f$
 */
double dsqrtdx(double x);

/*!\brief cosine function.
 *
 * The unit check is done with radians,
 * then the function double cos(double )
 * and double sin(double ) are used
 * to compute the values and uncertainty
 * values respectively.
 */
ParameterPhy cos(const ParameterPhy &rhs);

/*!\brief sine function.
 *
 * The unit check is done with radians,
 * then the function double sin(double )
 * and double cos(double ) are used
 * to compute the values and uncertainty
 * values respectively.
 */
ParameterPhy sin(const ParameterPhy &rhs);

/*!\brief arccosine function.
 *
 * The unit check is done with no unit,
 * then the function double acos(double )
 * and double dacosdx(double ) are used
 * to compute the values and uncertainty
 * values respectively.
 */
ParameterPhy acos(const ParameterPhy &rhs);

/*!\brief arcsine function.
 *
 * The unit check is done with no unit,
 * then the function double asin(double )
 * and double dasindx(double ) are used
 * to compute the values and uncertainty
 * values respectively.
 */
ParameterPhy asin(const ParameterPhy &rhs);

/*!\brief tangent function.
 *
 * The unit check is done with radians,
 * then the function double tan(double )
 * and double dtandx(double ) are used
 * to compute the values and uncertainty
 * values respectively.
 */
ParameterPhy tan(const ParameterPhy &rhs);

/*!\brief arctangent function.
 *
 * The unit check is done with no unit,
 * then the function double atan(double )
 * and double datandx(double ) are used
 * to compute the values and uncertainty
 * values respectively.
 */
ParameterPhy atan(const ParameterPhy &rhs);

/*!\brief exponential function.
 *
 * The unit check is done with no unit,
 * then the function double exp(double )
 * and double exp(double ) are used
 * to compute the values and uncertainty
 * values respectively.
 */
ParameterPhy exp(const ParameterPhy &rhs);

/*!\brief natural logarithm function.
 *
 * The unit check is done with no unit,
 * then the function double log(double )
 * and double dlogdx(double ) are used
 * to compute the values and uncertainty
 * values respectively.
 */
ParameterPhy log(const ParameterPhy &rhs);

/*!\brief base 10 logarithm function.
 *
 * The unit check is done with no unit,
 * then the function double log10(double )
 * and double dlog10dx(double ) are used
 * to compute the values and uncertainty
 * values respectively.
 */
ParameterPhy log10(const ParameterPhy &rhs);

/*!\brief absolute value function.
 *
 * The unit needs no management, it is copied, 
 * the absolute values are computed using the 
 * double abs(double ) function and a coefficient
 * of one, the uncertainty are simply copied.
 */
ParameterPhy fabs(const ParameterPhy &rhs);

/*!\brief Error function
 *
 * The unit check is done with no unit,
 * then the function double erf(double )
 * and double derfdx(double ) are used
 * to compute the values and uncertainty
 * values respectively.
 */
ParameterPhy erf(const ParameterPhy &z);

/*! \brief Root square function
 *
 * The unit check is done with a Units object
 */
ParameterPhy sqrt(const ParameterPhy &x);

/*!\brief power function with an integer.
 *
 * The power function with an integer uses the
 * operator Units Units::Units &operator*=(int)
 * to deal with the units. Then elevate all the
 * values to the power given, and combine the uncertainties.
 * \f[
 *      u_{\mathrm{pow(\alpha)}}^2 = \left(\alpha\nu^{\alpha-1}\right)^2 u^2
 * \f]
 * \f$\nu\f$ being the computed value of the parameter.
 *
 * As it needs another function than the double f(double ),
 * we do not use the calculateFunction(D)Values methods.
 */
ParameterPhy pow(const ParameterPhy &rhs, int p);

/*!\brief power function with a double.
 *
 * Unit management with double is not supported, and is
 * physically questionnable. Thus this method will only accept
 * unitless ParameterPhy.
 * Then the computations are done, the combined uncertainties are given
 * by the relation.
 * \f[
 *      u_{\mathrm{pow(X,\alpha)}}^2 = \left(\alpha\nu^{\alpha-1}\right)^2 u^2
 * \f]
 * \f$\nu\f$ being the computed value of the parameter.
 */
ParameterPhy pow(const ParameterPhy &rhs, double p);

/*!\brief power function with a double.
 *
 * The unit check is done with no unit.
 * Then the computations are done, the combined uncertainties are given
 * by the relation.
 * \f[
 *      u_{\mathrm{pow(\alpha,X)}}^2 = \left(\ln(\alpha) \alpha^X\right)^2 u^2
 * \f]
 * \f$\nu\f$ being the computed value of the parameter.
 */
ParameterPhy pow(double b, const ParameterPhy &power);

/*!\brief power function with a ParameterPhy.
 *
 * Unit management with double is not supported, and is
 * physically questionnable. Thus this method will only accept
 * unitless ParameterPhy for the base, and as it should be in any
 * case, for the power.
 * Thus unit check is done both for the power and base with
 * no unit.
 * Then the computations are done, the combined uncertainties are given
 * by the relation.
 * \f[
 *      u_{\varphi_1^{\varphi_2}}^2 = 
 *              \left(\nu_2\nu_1^{\nu_2-1}\right)^2 u_1^2 +
 *              \left(\ln(\nu_1)\nu_1^{\nu_2}\right)^2 u_2^2
 * \f]
 * \f$\nu_i\f$ being the computed value of the parameter \f$i\f$.
 *
 * A complementary check must be performed before the computations:
 * the base and power must have the same number of values, or at least
 * one of them should have only one value. If not, it sends an error.
 */
ParameterPhy pow(const ParameterPhy &base, const ParameterPhy &p);
}

#endif
