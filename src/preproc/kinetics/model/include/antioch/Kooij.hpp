#ifndef _KOOIJ_KINETICS_MODELS_
#define _KOOIJ_KINETICS_MODELS_

#include "antioch/HercourtEssen.hpp"
#include "antioch/Arrhenius.hpp"

namespace Antioch{
 /*!\file Kooij.hpp
 * \brief contains the corresponding class
 *
 * \class Kooij
 * \ingroup kinmod
 * \brief The Kooij kinetics model
 *
 * The Kooij equation is
 * \f[
 *      \KooijEquation
 * \f]
 * There is a biais in the Kooij uncertainty calculation:
 * if we decompose the uncertainty calculation, as is done in the program, we obtain, for the Kooij
 * equation:
 * \f[
 * \renewcommand{\arraystretch}{1.5}
 * \begin{array}{r@{\;=\;}ll}\toprule
 *      u_{\Kooij}^2 & u_{g}^2 f^2 + u_{f}^2 g^2 & \mathrm{with}\quad
 *                                                                      \begin{array}{l}
 *                                                                       f = A \left(\frac{T}{\To}\right)^\beta \\
 *                                                                       g = \exp\left(-\frac{E_a}{\R T}\right) 
 *                                                                       \end{array}\tabularnewline\cmidrule(rl){1-3}\addlinespace
 *      u_{f}^2  & p^2 u_A^2 + u_{p}^2 A^2            & \mathrm{with}\quad
 *                                                                      p = \left(\frac{T}{\To}\right)^\beta \tabularnewline
 *      u_{p}^2  & \left(\beta q^{\beta-1}\right)^2 u_{q}^2 + \left(\ln(q) q^\beta\right)^2 u_\beta^2
 *                                                             & \mathrm{with}\quad
 *                                                                      q = \frac{T}{\To} \tabularnewline
 *      u_{q}^2  & \left(\frac{1}{\To}\right)^2 u_T^2 + \left(-\frac{T}{\To^2}\right)^2 u_{\To}^2 \tabularnewline
 *
 *      \addlinespace[10pt]\cmidrule{2-2}\addlinespace
 *      
 *      u_{g}^2  & \left(\exp\left(\varphi\right)\right)^2 u_{\varphi}^2
 *                                                             & \mathrm{with}\quad
 *                                                                      \varphi = -\frac{E_a}{\R T} \tabularnewline
 *      u_{\varphi}^2
 *                  & \left(\frac{E_a}{v^2}\right)^2 u_{v}^2 + \left(-\frac{1}{v}\right)^2 u_{E_a}^2
 *                                                             & \mathrm{with}\quad
 *                                                                      v = \R T \tabularnewline
 *      u_{v}^2  & \R^2 u_T^2 + T^2 u_{\R}^2
 *  \tabularnewline
 *  \bottomrule                                                                   
 *  \end{array}
 * \f]
 * Finally, decomposed in each of the variables, we obtain:
 * \f[
 * \renewcommand{\arraystretch}{1.5}
 * \begin{array}{rll}
 * u_{\Kooij}^2 =& \left(\KooijEq\right)^2 \times 
 * \tabularnewline
 *     ( \tabularnewline
 *      &    \left[\left(\frac{E_a}{\R T}\right)^2 + \beta^2\right]\frac{1}{T^2}   & u_T^2     \tabularnewline
 *     &+    \left(\frac{E_a}{\R^2 T}\right)^2                                     & u_{\R}^2  \tabularnewline
 *     &+    \left(-\frac{1}{\R T}\right)^2                                        & u_{E_a}^2 \tabularnewline
 *     &+    \left(\frac{1}{A}\right)^2                                            & u_A^2     \tabularnewline
 *     &+    \left(\beta\frac{1}{\To}\right)^2                                     & u_{\To}^2 \tabularnewline
 *     &+    \left(\ln\left(\frac{T}{\To}\right)\right)^2                          & u_\beta^2 \tabularnewline
 * )
 * \end{array}
 * \f]
 * If we apply the GUM formulae \f$\left(\displaystyle\sum_i\left(\frac{\partial f}{\partial x_i}\right)^2u_{x_i}^2\right)\f$,
 * we obtain the decomposition:
 * \f[
 * \renewcommand{\arraystretch}{1.5}
 * \begin{array}{rll}
 * u_{\Kooij}^2 =& \left(\KooijEq\right)^2 \times 
 *    \tabularnewline ( \tabularnewline
 *               &  \left(\frac{1}{A}\right)^2                                                                     & u_A^2     \tabularnewline
 *               &+ \left[\left(\frac{E_a}{\R T}\right)^2 + \beta^2 + 2\beta\frac{E_a}{\R T}\right]\frac{1}{T^2}   & u_T^2     \tabularnewline
 *               &+ \left(\frac{E_a}{\R^2 T}\right)^2                                                              & u_{\R}^2  \tabularnewline
 *               &+ \left(-\frac{1}{\R T}\right)^2                                                                 & u_{E_a}^2 \tabularnewline
 *               &+ \left(\beta\frac{1}{\To}\right)^2                                                              & u_{\To}^2 \tabularnewline
 *               &+ \left(\ln\left(\frac{T}{\To}\right)\right)^2                                                   & u_\beta^2 \tabularnewline
 * )
 * \end{array}
 * \f]
 * Thus, the biais over the variance on a Kooij equation is
 * \f[
 *      \text{biais}_{\Kooij} = \KooijBiais
 * \f]
 *
 * The default value for \f$\To\f$ is 300 K. The temperature at which
 * the rate constant should be calculated can be given in two ways:
 *   - by setting the ParameterPhy pointer T to the wanted temperature,
 *   - by giving the temperature in a double form. If the parameters
 *      contains more than one value, the index to the wanted value
 *      should be given. By default, the calculations are performed
 *      with the first value.
 *
 * As it inherits multiple times KineticsModel, the preexponential factor is
 * stored somewhere else.
 */
class Kooij:public HercourtEssen,
            public Arrhenius
{
   public:
/*!\brief Default constructor*/
     Kooij(){}
/*!\brief Copy constructor*/
     Kooij(const Kooij &rhs) {*this = rhs;}
/*!\brief Constructor*/
     Kooij(const ParameterPhy &Pre, const ParameterPhy &pow, const ParameterPhy &aEn, const ParameterPhy &Tref = PhyCon::TrefKin):
        HercourtEssen(Pre,pow,Tref),
        Arrhenius(ParameterPhy("one",1.0,0.0,CORE_UNCERTAINTY_TYPE_NONE,""),aEn){}
/*!\brief Automatic choosing constructor*/
     Kooij(const std::vector<ParameterPhy> &pars);
/*!\brief Destructor*/
     ~Kooij(){}

/*!\brief Initialize from a vector of parameters*/
     void init(const std::vector<ParameterPhy> &pars);
/*!\brief Change the values from a vector of double*/
     void reset(const std::vector<double> &pars);
/*!\brief Rate constant calculator
 *
 * It checks first if the parameters (\f$A\f$, \f$\beta\f$ and \f$E_a\f$)
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
     double getRateConstantT(double Tem, int i = 0) const;//check parameters ok

/*!\brief Default configuration: unit*/
     const std::string getDefaultUnit(int i)       const;
/*!\brief Default configuration: min*/
     double getDefaultMin(int i)                   const;
/*!\brief Default configuration: max*/
     double getDefaultMax(int i)                   const;
/*!\brief Default configuration: all in ParameterPhy*/
     const ParameterPhy getDefaultParameter(int i) const;

/*!\brief Parameters getter*/
     const std::vector<ParameterPhy> getParameters() const;
/*!\brief Number of parameters getter*/
     unsigned int nPars()       const {return 3;}
/*!\brief Number of optionnal parameters getter*/
     unsigned int nOptionPars() const {return 1;}


/*!\brief Assignement operator*/
     Kooij &operator=(const Kooij &rhs);

   private:

/*!\brief showKinModel()*/
     void showKinModel(std::ostream &out = std::cout) const;

/*\brief Check all parameters*/
     void checkParams() const;
/*\brief Sort a vector of parameters
 *
 * The parameters are sorted in this order:
 *  - the parameter having K as unit or being homogeneous to J/mol is the
 *      activation energy;
 *  - the parameter found with no unit is \f$\beta\f$;
 *  - if found, the parameter found with exactly K unit is \f$\To\f$,
 *    if not found, \f$\To\f$ is set to default (300 K);
 *  - the remaining parameter is the pre-exponential factor.
 */
     std::vector<const ParameterPhy*> whatParametersAreWe(const std::vector<ParameterPhy> &pars);
/*\brief Templated equation
 *
 * This is the only place where the equation is written,
 * used in the ParameterPhy and double version.
 */
     template<class M>
     M kEquation(const M& temp, const M& pre, const M& beta, const M& ea, const M& tref) const;
};
}
#endif
