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
#ifndef _ARRHENIUS_KINETICS_MODELS_
#define _ARRHENIUS_KINETICS_MODELS_

#include "antioch/KineticsModel.hpp"

namespace Antioch{
 /*!\file Arrhenius.hpp
 * \brief contains the corresponding class
 *
 * \class Arrhenius
 * \ingroup kinmod
 * \brief The Arrhenius kinetics model
 *
 * The Arrhenius equation is:
 * \f[
 *     \ArrhEquation
 * \f]
 * 
 * The uncertainty following the GUM recommandation is:
 * \f[
 * \renewcommand{\arraystretch}{1.5}
 * \begin{array}{rll}
 *    u_{\Arrh}^2  =& \left(\ArrhEq\right)^2\times \tabularnewline
 *    ( \tabularnewline
 *                  &   \left(\frac{1}{A}\right)^2        & u_{A}^2   \tabularnewline
 *                  &+  \left(-\frac{1}{\R T}\right)^2    & u_{E_a}^2 \tabularnewline
 *                  &+  \left(\frac{E_a}{\R T^2}\right)^2 & u_{T}^2   \tabularnewline
 *                  &+  \left(\frac{E_a}{T \R^2}\right)^2 & u_{\R}^2  \tabularnewline
 *    )
 * \end{array}
 * \f]
 * 
 * There is no biais in the uncertainty calculation.
 *
 * The temperature at which
 * the rate constant should be calculated can be given in two ways:
 *   - by setting the ParameterPhy pointer T to the wanted temperature,
 *   - by giving the temperature in a double form. If the parameters
 *      contains more than one value, the index to the wanted value
 *      should be given. By default, the calculations are performed
 *      with the first value.
 *
 *  This class use the reduced exponential parameter (\f$E_a / \R\f$) instead
 *  of the activation energy.
 */
class Arrhenius:virtual public KineticsModel
{
   public:
/*!\brief Default constructor*/
     Arrhenius(const ParameterPhy &ea = ParameterPhy("Activation energy","J/mol")):Ea(ea){}
/*!\brief Constructor*/
     Arrhenius(const ParameterPhy &Pre, const ParameterPhy &aEn):
        KineticsModel(Pre),Ea(aEn){setExpPa();}
/*!\brief Copy constructor*/
     Arrhenius(const Arrhenius &rhs) {*this = rhs;}
/*!\brief Automatic choosing constructor*/
     Arrhenius(const std::vector<ParameterPhy> &pars);
/*!\brief Destructor*/
    virtual ~Arrhenius(){}

/*!\brief Initialize from a vector of parameters*/
     virtual void init(const std::vector<ParameterPhy> &pars);
/*!\brief Change the values from a vector of double*/
     virtual void reset(const std::vector<double> &pars);
/*!\brief Rate constant calculator
 *
 * It checks first if the parameters (\f$A\f$ and \f$E_a\f$)
 * contains values, if the temperature is defined, and returns the
 * ParameterPhy associated with the rate constant.
 */
     virtual const ParameterPhy rateConstant(int nk = 0);//check parameters and pointer ok
/*!\brief Rate constant calculator
 *
 * It checks if the parameters are defined, and
 * calculate the rate constant, using value \f$i\f$
 * of the parameter (0 by default), at temperature Tem.
 */
     virtual double getRateConstantT(double Tem, int i = 0) const;//check parameters ok

/*!\brief \f$E_a\f$ getter*/
     ParameterPhy getEa()   const {return Ea;}

/*!\brief \f$E_a\f$ setter*/
     void setEa(const ParameterPhy &ae);
/*!\brief \f$E_a\f$ setter*/
     void resetEa(double ae);

/*!\brief Parameters getter*/
     const std::vector<ParameterPhy> getParameters() const;
/*!\brief Number of parameters getter*/
     virtual unsigned int nPars() const {return 2;}

/*!\brief Default configuration: unit*/
     const std::string getDefaultEaUnit() const {return "J/mol";}
/*!\brief Default configuration: min*/
     double getDefaultEaMin()             const {return -5e4;}
/*!\brief Default configuration: max*/
     double getDefaultEaMax()             const {return 5e4;}
/*!\brief Default configuration: all in ParameterPhy*/
     const ParameterPhy getDefaultEa()    const;

/*!\brief Default configuration: unit*/
     virtual const std::string getDefaultUnit(int i)       const;
/*!\brief Default configuration: min*/
     virtual double getDefaultMin(int i)                   const;
/*!\brief Default configuration: max*/
     virtual double getDefaultMax(int i)                   const;
/*!\brief Default configuration: all in ParameterPhy*/
     virtual const ParameterPhy getDefaultParameter(int i) const;

/*!\brief Assignement operator*/
     Arrhenius &operator=(const Arrhenius &rhs);

   protected:
/*\brief \f$E_a\f$ parameter*/
     ParameterPhy Ea;
/*\brief Reduced exponential parameter*/
     ParameterPhy expPa;
/*\brief Reduced exponential parameter used for calculations*/
     ParameterPhy expPaCal;

/*\brief Check all parameters*/
     virtual void checkParams() const;
/*\brief Sort a vector of parameters
 *
 * The sorting is done unit-wise.
 * First it looks for Ea, units accepted are:
 *    - homogeneous to J/mol (actual activation energy),
 *    - K (reduced activation energy)
 * then whatever remains is the preexponential factor, no
 * assumption on its unit is made.
 */
     virtual std::vector<const ParameterPhy*> whatParametersAreWe(const std::vector<ParameterPhy> &pars);

   private:
/*\brief Sets the reduced exponential parameter*/
     void setExpPa();
/*\brief Templated equation
 *
 * This is the only place where the equation is written,
 * used in the ParameterPhy and double version.
 */
     template<class M>
     M kEquation(const M& temp, const M& pre, const M& ea) const;
};
}
#endif
