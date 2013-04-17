//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
#ifndef _BERTHELOT_HERCOURT_ESSEN_KINETICS_MODELS_
#define _BERTHELOT_HERCOURT_ESSEN_KINETICS_MODELS_

#include "antioch/HercourtEssen.hpp"
#include "antioch/Berthelot.hpp"

namespace Antioch{
 /*!\file BerthelotHercourtEssen.hpp
 * \brief contains the corresponding class
 *
 * \class BerthelotHercourtEssen
 * \ingroup kinmod
 * \brief The Berthelot Hercourt-Essen kinetics model
 *
 * The Berthelot Hercourt-Essen equation is:
 * \f[
 *     \BertHerEssEquation
 * \f]
 * There is a biais in the calculation. Following the program protocol, we have:
 * \f[
 * \renewcommand{\arraystretch}{1.5}
 * \begin{array}{r@{\;=\;}ll}\toprule
 *      u_{\BertHerEss}^2 & u_{g}^2 f^2 + u_{f}^2 g^2 & \mathrm{with}\quad
 *                                                                      \begin{array}{l}
 *                                                                       f = A \left(\frac{T}{\To}\right)^\beta \\
 *                                                                       g = \exp\left(DT\right) 
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
 *                                                                      \varphi = DT \tabularnewline
 *      u_{\varphi}^2
 *                  & D^2 u_{T}^2 + T^2 u_{D}^2
 *  \tabularnewline
 *  \bottomrule                                                                   
 *  \end{array}
 * \f]
 * which gives for each variable:
 * \f[
 * \renewcommand{\arraystretch}{1.5}
 * \begin{array}{rll}
 * u_{\BertHerEss}^2 =& \left(\BertHerEssEq\right)^2 \times 
 *    \tabularnewline ( \tabularnewline
 *                    &  \left(\frac{1}{A}\right)^2                        & u_A^2     \tabularnewline
 *                    &+ \left[D^2 + \left(\frac{\beta}{T}\right)^2\right] & u_T^2     \tabularnewline
 *                    &+ T^2                                               & u_{D}^2   \tabularnewline
 *                    &+ \left(\frac{\beta}{\To}\right)^2                  & u_{\To}^2 \tabularnewline
 *                    &+ \left(\ln\left(\frac{T}{\To}\right)\right)^2      & u_\beta^2 \tabularnewline
 * )
 * \end{array}
 * \f]
 * According to the GUM, we obtain:
 * \f[
 * \renewcommand{\arraystretch}{1.5}
 * \begin{array}{rll}
 * u_{\BertHerEss}^2 =& \left(\BertHerEssEq\right)^2 \times 
 *    \tabularnewline ( \tabularnewline
 *                    &  \left(\frac{1}{A}\right)^2                                            & u_A^2     \tabularnewline
 *                    &+ \left[D^2 + \left(\frac{\beta}{T}\right)^2 + 2D\frac{\beta}{T}\right] & u_T^2     \tabularnewline
 *                    &+ T^2                                                                   & u_{D}^2   \tabularnewline
 *                    &+ \left(\frac{\beta}{\To}\right)^2                                      & u_{\To}^2 \tabularnewline
 *                    &+ \left(\ln\left(\frac{T}{\To}\right)\right)^2                          & u_\beta^2 \tabularnewline
 * )
 * \end{array}
 * \f]
 * thus, the biais is:
 * \f[
 *    \text{biais}_{\BertHerEss} = \BertHerEssBiais
 * \f]
 *
 * The temperature at which
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
class BerthelotHercourtEssen:public Berthelot,
                             public HercourtEssen
{
   public:
/*!\brief Default constructor*/
     BerthelotHercourtEssen(){}
/*!\brief Constructor*/
     BerthelotHercourtEssen(const ParameterPhy &Pre, const ParameterPhy &pow, const ParameterPhy &d, const ParameterPhy &Tref = PhyCon::TrefKin):
        Berthelot(Pre,d),
        HercourtEssen(Pre,pow,T0){}
/*!\brief Automatic choosing constructor*/
     BerthelotHercourtEssen(const std::vector<ParameterPhy> &pars);
/*!\brief Destructor*/
     ~BerthelotHercourtEssen(){}

/*!\brief Initialize from a vector of parameters*/
     void init(const std::vector<ParameterPhy> &pars);
/*!\brief Change the values from a vector of double*/
     void reset(const std::vector<double> &pars);
/*!\brief Rate constant calculator
 *
 * It checks first if the parameters (\f$A\f$, \f$\beta\f$ and \f$D\f$)
 * contains values, if the temperature is defined, and returns the
 * ParameterPhy associated with the rate constant.
 */
     const ParameterPhy rateConstant(int kn = 0);//check parameters and pointer ok
/*!\brief Rate constant calculator
 *
 * It checks if the parameters are defined, and
 * calculate the rate constant, using value \f$i\f$
 * of the parameter (0 by default), at temperature Tem.
 */
     double getRateConstantT(double Tem, int i = 0) const;

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
     unsigned int nPars()       const {return 2;}
/*!\brief Number of optionnal parameters getter*/
     unsigned int nOptionPars() const {return 1;}

/*!\brief Assignement operator*/
     BerthelotHercourtEssen &operator=(const BerthelotHercourtEssen &rhs);

   private:

/*\brief Check all parameters*/
     void checkParams() const;
/*\brief Sort a vector of parameters
 *
 * The parameters are sorted in this order:
 *  - the parameter having K-1 as unit is the exponential parameter;
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
     M kEquation(const M& temp , const M& pre, const M& beta, const M& d, const M& Tref) const;
};
}
#endif
