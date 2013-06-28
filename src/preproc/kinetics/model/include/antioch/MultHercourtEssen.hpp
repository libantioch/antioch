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
#ifndef _MULTIPLE_HERCOURT_ESSEN_KINETICS_MODELS_
#define _MULTIPLE_HERCOURT_ESSEN_KINETICS_MODELS_

#include "antioch/KineticsModel.hpp"
#include "antioch/physical_constants.hpp"

namespace Antioch{
/*!\file MultHercourtEssen.hpp
 * \brief contains the corresponding class
 *
 * \class MultHercourtEssen
 * \ingroup kinmod
 * \brief Kinetics model for multiple Hercourt-Essen equation
 *
 * The multiple Hercourt-Essen equation is:
 * \f[
 *    \mHerEssEquation
 * \f]
 *
 * The uncertainty following the GUM recommandation is:
 * \f[
 * \renewcommand{\arraystretch}{1.5}
 * \begin{array}{rll}
 *    u_{\mHerEss}^2 =& \left(\mHerEssEq\right)^2\times \tabularnewline
 *    ( \tabularnewline
 *                  &   \left(\frac{1}{A}\right)^2                       & u_{A}^2   \tabularnewline
 *                  &+  \left(\ln\left(\frac{T_1}{\Too}\right)\right)^2  & u_{\beta_1}^2 \tabularnewline
 *                  &+  \left(\ln\left(\frac{T_2}{\Tto}\right)\right)^2  & u_{\beta_2}^2 \tabularnewline
 *                  &+  \left(\frac{\beta_1}{T_1}\right)^2               & u_{T_1}^2 \tabularnewline
 *                  &+  \left(\frac{\beta_2}{T_2}\right)^2               & u_{T_2}^2 \tabularnewline
 *                  &+  \left(-\beta_1 \Too\right)^2                     & u_{\Too}^2 \tabularnewline
 *                  &+  \left(-\beta_2 \Tto\right)^2                     & u_{\Tto}^2 \tabularnewline
 *    )
 * \end{array}
 * \f]
 * 
 * There is no biais in the uncertainty calculation, we consider that \f$\Too\f$ and \f$\Tto\f$ are
 * two different parameters.
 *
 * The default value for \f$\Too\f$ and \f$\Tto\f$ is 300 K. The temperature at which
 * the rate constant should be calculated can be given in two ways:
 *   - by setting the ParameterPhy pointers T and Ti to the wanted temperature,
 *   - by giving the temperatures in a double form. If the parameters
 *      contains more than one value, the index to the wanted value
 *      should be given. By default, the calculations are performed
 *      with the first value.
 */
class MultHercourtEssen:public KineticsModel
{
   public:
/*!\brief Default constructor*/
     MultHercourtEssen(const ParameterPhy &be = ParameterPhy("Power parameter 1",""),
                       const ParameterPhy &bi = ParameterPhy("Power parameter 2","")):
                        betae(be),betaeCal(be),betai(bi),betaiCal(bi),
                        T0e(PhyCon::TrefKin),T0i(PhyCon::TrefKin){}
/*!\brief Copy constructor*/
     MultHercourtEssen(const MultHercourtEssen &rhs) {*this = rhs;}
/*!\brief Constructor*/
     MultHercourtEssen(const ParameterPhy &Pre, const ParameterPhy &powe, const ParameterPhy &powi, 
        const ParameterPhy &Trefe = PhyCon::TrefKin, const ParameterPhy &Trefi = PhyCon::TrefKin):
        KineticsModel(Pre),betae(powe),betaeCal(powe),betai(powi),betaiCal(powi),T0e(Trefe),T0i(Trefi),
        Ti(NULL){}
/*!\brief Automatic choosing constructor*/
     MultHercourtEssen(const std::vector<ParameterPhy> &pars);
/*!\brief Destructor*/
     ~MultHercourtEssen(){}

/*!\brief Initialize from a vector of parameters*/
     void init(const std::vector<ParameterPhy> &pars);
/*!\brief Change the values from a vector of double*/
     void reset(const std::vector<double> &pars);
/*!\brief Rate constant calculator
 *
 * It checks first if the parameters (\f$A\f$, \f$\beta_1\f$ and \f$\beta_2\f$)
 * contains values, if the temperature is defined, and returns the
 * ParameterPhy associated with the rate constant.
 */
     const ParameterPhy rateConstant(int nk = 0);//check parameters and pointer ok
/*!\brief Rate constant calculator
 *
 * It checks if the parameters are defined, and
 * calculate the rate constant, using value \f$i\f$
 * of the parameter (0 by default), at temperature Tem.
 */
     double getRateConstantT(double Tem, double Tim, int i = 0) const;//check parameters ok
/*!\brief Unpurify this class*/
     double getRateConstantT(double Tem, int i = 0) const;//sends back an error message

/*!\brief \f$\beta_1\f$ getter*/
     ParameterPhy getBetae() const {return betae;}
/*!\brief \f$\beta_2\f$ getter*/
     ParameterPhy getBetai() const {return betai;}
/*!\brief \f$\Too\f$ getter*/
     ParameterPhy getT0e()   const {return T0e;}
/*!\brief \f$\Tto\f$ getter*/
     ParameterPhy getT0i()   const {return T0i;}

/*!\brief \f$T\f$ setter*/
     void setT(ParameterPhy *Tused, int nT = 0);
/*!\brief \f$T\f$ getter*/
     ParameterPhy * getT(int nT = 0) const;

/*!\brief \f$\beta_1\f$ setter*/
     void setBetae(const ParameterPhy &pow) {betae.replace(pow);betaeCal.replace(pow);}
/*!\brief \f$\beta_2\f$ setter*/
     void setBetai(const ParameterPhy &pow) {betai.replace(pow);betaiCal.replace(pow);}
/*!\brief \f$\beta_1\f$ setter*/
     void resetBetae(double pow)            {resetPar(betae,pow);resetPar(betaeCal,pow);}
/*!\brief \f$\beta_2\f$ setter*/
     void resetBetai(double pow)            {resetPar(betai,pow);resetPar(betaiCal,pow);}
/*!\brief \f$\Too\f$ setter*/
     void setT0e(const ParameterPhy &Tref)  {T0e.replace(Tref);}
/*!\brief \f$\Tto\f$ setter*/
     void setT0i(const ParameterPhy &Tref)  {T0i.replace(Tref);}

/*!\brief Default configuration: unit*/
     const std::string getDefaultBetaeUnit() const {return std::string();}
/*!\brief Default configuration: min*/
     double getDefaultBetaeMin() const {return 0.;}
/*!\brief Default configuration: max*/
     double getDefaultBetaeMax() const {return 2.;}
/*!\brief Default configuration: all in ParameterPhy*/
     const ParameterPhy getDefaultBetae() const;
/*!\brief Default configuration: unit*/
     const std::string getDefaultBetaiUnit() const {return std::string();}
/*!\brief Default configuration: min*/
     double getDefaultBetaiMin() const {return 0.;}
/*!\brief Default configuration: max*/
     double getDefaultBetaiMax() const {return 2.;}
/*!\brief Default configuration: all in ParameterPhy*/
     const ParameterPhy getDefaultBetai() const;

/*!\brief Default configuration: unit*/
     virtual const std::string getDefaultUnit(int i) const;
/*!\brief Default configuration: min*/
     virtual double getDefaultMin(int i) const;
/*!\brief Default configuration: max*/
     virtual double getDefaultMax(int i) const;
/*!\brief Default configuration: all in ParameterPhy*/
     virtual const ParameterPhy getDefaultParameter(int i);



/*!\brief First temperature pointer setter*/
     void setTe(ParameterPhy *Teused)  {T = Teused;}
/*!\brief First temperature getter*/
     ParameterPhy * getTe()      const {return T;}
/*!\brief Second temperature pointer setter*/
     void setTi(ParameterPhy *Tiused)  {Ti = Tiused;}
/*!\brief Second temperature getter*/
     ParameterPhy * getTi()      const {return Ti;}
/*!\brief Test to see if Ti is defined*/
     bool isTiDefined()          const {return !(Ti == NULL);}

/*!\brief Parameters getter*/
     const std::vector<ParameterPhy> getParameters() const;
/*!\brief Number of parameters getter*/
     unsigned int nPars()       const {return 3;}
/*!\brief Number of optionnal parameters getter*/
     unsigned int nOptionPars() const {return 2;}

/*!\brief Assignement operator*/
     MultHercourtEssen &operator=(const MultHercourtEssen &rhs);

   private:
/*\brief \f$\beta_1\f$ parameter*/
     ParameterPhy betae;
/*\brief \f$\beta_1\f$ parameter used for calculations*/
     ParameterPhy betaeCal;
/*\brief \f$\beta_2\f$ parameter*/
     ParameterPhy betai;
/*\brief \f$\beta_2\f$ parameter used for calculations*/
     ParameterPhy betaiCal;
/*\brief \f$\Too\f$ parameter*/
     ParameterPhy T0e;
/*\brief \f$\Tto\f$ parameter*/
     ParameterPhy T0i;

/*\brief pointer to second temperature*/
     ParameterPhy *Ti;

/*\brief Check all parameters*/
     void checkParams() const;
/*\brief Sort a vector of parameters
 *
 * The parameters are sorted in this order:
 *  - the first parameter found with no unit is \f$\beta_1\f$;
 *  - the second parameter found with no unit is \f$\beta_2\f$;
 *  - if found, the first parameter found with exactly K unit is \f$\Too\f$,
 *    if not found, \f$\Too\f$ and \f$\Tto\f$ are set to default (300 K);
 *  - if found, the second parameter found with exactly K unit is \f$\Tto\f$,
 *    if not found and, \f$\Tto\f$ is set to \f$\Too\f$.
 *  - the remaining parameter is the pre-exponential factor.
 */
     std::vector<const ParameterPhy*> whatParametersAreWe(const std::vector<ParameterPhy> &pars);
/*\brief Templated equation
 *
 * This is the only place where the equation is written,
 * used in the ParameterPhy and double version.
 */
     template<class M>
     M kEquation(const M& temp1, const M& temp2, const M& pre, const M& beta1, const M& beta2, const M& Tref1, const M& Tref2) const;
};
}
#endif
