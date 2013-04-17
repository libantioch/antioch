//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
#ifndef _FALLOFF_REACTION_
#define _FALLOFF_REACTION_

#include "antioch/Process.hpp"

namespace Antioch{
/*!\file FalloffProcess.hpp
 * \brief Falloff base class
 * \todo finish the documentation, the biais calculations
 *
 * \class FalloffProcess
 * \brief base class for falloff processes
 *
 * A falloff process is a reaction having two regimes,
 * a high pressure (\f$k_\infty(T)\f$ and a low pressure 
 * (\f$k_0(T) \ce{[M]}\f$) regime. The falloff expression is the
 * mathematical function that, using those two regimes, will
 * describe the rate constant at any pressure, i.e., 
 * \f$k(T,P) = f(k_0(T)\ce{[M]},k_\infty(T))\f$.
 * The general expression is given by:
 * \f[
 *  k(T,P) = k_\infty \frac{P_r}{1 + P_r} F
 * \f]
 * with
 * \f[
 *   P_r = \frac{\ce{[M]} k_0}{k_\infty}
 * \f]
 * with \f$\ce{[M]}\f$ the concentration of the mixture, being thus
 * the pressure indicator. As \f$P_r\f$
 * should be unitless, the unit of \f$k_0\ce{[M]}\f$ is the same as
 * \f$k_\infty\f$, thus \f$k_0\f$ is a rate constant of an order of
 * 1 less than \f$k_\infty\f$.
 *
 * Thus if \f$k_\infty\f$ is a rate constant of order 2, its unit
 * is \f$\text{cm}^3\text{s}^{-1}\f$, thus \f$k_0\f$ will have a
 * unit of \f$\text{s}^{-1}\f$
 *
 * For this class, the rate constant provided by the class Process
 * is considered the low pressure \f$k_0\f$ constant.
 *
 * There is a biais in the rate constant uncertainty calculation.
 *
 * The correction on the variance is:
 * \f[
 *   \text{biais} = \left(-2  \left(k_\infty F\right)^2 \frac{P_r}{\left(1 + P_r\right)^3}\right)^2 u^2_{P_r}
 * \f]
 *
 * This class is ab abstract class so every falloff model has to
 * implement the \f$F\f$ function, in ParameterPhy and double version
 * at least.
 */

class FalloffProcess:public Process
{
  public:
/*!\brief Default constructor*/
      FalloffProcess(){rates.push_back(NULL);}
/*!\brief Full constructor*/
      FalloffProcess(const std::string &kinMod, const std::vector<ParameterPhy> &pars0, const std::vector<ParameterPhy> &parsinf, ParameterPhy *m = NULL);
/*!\brief Destructor*/
      ~FalloffProcess();

/*!\brief Initialization
 *
 * 1 - set kinetics model
 * 2 - set parameters k0
 * 3 - set parameters kinf
 */
        virtual void init(const std::string &kinMod, const std::vector<ParameterPhy> &pars0, const std::vector<ParameterPhy> &parsinf, ParameterPhy *m = NULL);
/*!\brief Sets the temperature*/
        virtual void setTemperature(ParameterPhy *t);
/*!\brief Sets F parameters*/
        virtual void setFalloffParameters(const std::vector<ParameterPhy> &pars){} //nothing to do here
/*!\brief Resets F parameters*/
        virtual void resetParameters(const std::vector<double> &pars);

/*!\brief All parameters getter*/
        virtual const std::vector<ParameterPhy> getParameters() const;

/*!\brief Rate constant calculator
 *
 * Temperature is given, the mixture concentration is defined:
 *   - M
 */
        virtual ParameterPhy rateConstantT(double t, int i = 0);
/*!\brief Rate constant calculator
 *
 * Mixture concentration is given, the temperature is defined:
 *   - T
 */
        virtual ParameterPhy rateConstantM(double m, int i = 0);
/*!\brief Rate constant calculator
 *
 * All is provided.
 */
        virtual double getRateConstantTM(double t, double m, int i = 0) const;

  protected:
/*!\brief Virtual pure function for double version of F
 *
 * Here M is defined, T is provided
 */
     virtual ParameterPhy logFT(double t, int i = 0) = 0;
/*!\brief Virtual pure function for double version of F
 *
 * Here T is defined, M is provided
 */
     virtual ParameterPhy logFM(double m, int i = 0) = 0;
/*!\brief Virtual pure function for double version of F
 *
 * Here all (T and M) is provided
 */
     virtual double logFTM(double t, double m, int i = 0) const = 0;

/*\brief reset the parameters for F*/
     virtual void resetF(std::vector<double> &pars) = 0;

/*!\brief Reduced pressure, double version
 *
 * T and M are provided.
 */
     double PrTM(double t, double m, int i = 0) const;
/*!\brief Reduced pressure, double version
 *
 * T is provided, M is defined.
 */
     ParameterPhy PrT(double t, int i = 0) const;
/*!\brief Reduced pressure, double version
 *
 * M is provided, T is defined.
 */
     ParameterPhy PrM(double m, int i = 0) const;

  private:
     template<class N>
     N PrEquation(const N &lowpres, const N &highpres, const N& conc) const;
     template<class N>
     N kEquation(const N &highpres, const N &pr, const N& f) const;
};
}
#endif
