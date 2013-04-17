//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
#ifndef _HERCOURT_ESSEN_KINETICS_MODELS_
#define _HERCOURT_ESSEN_KINETICS_MODELS_

#include "antioch/KineticsModel.hpp"
#include "antioch/physical_constants.hpp"
namespace Antioch{

/*!\file HercourtEssen.hpp
 * \brief Contains the corresponding class
 *
 * \class HercourtEssen
 * \ingroup kinmod
 * \brief Kinetics model for Hercourt-Essen equation
 *
 *  The Hercourt-Essen equation is:
 * \f[
 *     \HerEssEquation
 * \f]
 * The uncertainty calculated following the GUM recommandation is:
 * \f[
 * \renewcommand{\arraystretch}{1.5}
 * \begin{array}{rll}
 *    u_{\HerEss}^2 =& \left(\HerEssEq\right)^2\times \tabularnewline
 *    ( \tabularnewline
 *                   &   \left(\frac{1}{A}\right)^2                   & u_{A}^2     \tabularnewline
 *                   &+  \left(\ln\left(\frac{T}{\To}\right)\right)^2 & u_{\beta}^2 \tabularnewline
 *                   &+  \left(\frac{\beta}{T}\right)^2               & u_{T}^2     \tabularnewline
 *                   &+  \left(-\frac{\beta}{\To}\right)^2            & u_{\To}^2   \tabularnewline
 *    )
 * \end{array}
 * \f]
 *
 * There is no uncertainty biais for this kinetics model.
 *
 * The default value for \f$\To\f$ is 300 K. The temperature at which
 * the rate constant should be calculated can be given in two ways:
 *   - by setting the ParameterPhy pointer T to the wanted temperature,
 *   - by giving the temperature in a double form. If the parameters
 *      contains more than one value, the index to the wanted value
 *      should be given. By default, the calculations are performed
 *      with the first value.
 */
class HercourtEssen:virtual public KineticsModel
{
   public:
/*!\brief Default constructor*/
     HercourtEssen(const ParameterPhy &pow = ParameterPhy("Power parameter","")):beta(pow),betaCal(pow),T0(PhyCon::TrefKin){}
/*!\brief Constructor*/
     HercourtEssen(const ParameterPhy &Pre, const ParameterPhy &pow, const ParameterPhy &Tref = PhyCon::TrefKin):
        KineticsModel(Pre),beta(pow),betaCal(pow),T0(Tref){}
/*!\brief Copy constructor*/
     HercourtEssen(const HercourtEssen &rhs) {*this = rhs;}
/*!\brief Automatic choosing constructor*/
     HercourtEssen(const std::vector<ParameterPhy> &pars);
/*!\brief Destructor*/
    virtual  ~HercourtEssen(){}

/*!\brief Initialize from a vector of parameters*/
     virtual void init(const std::vector<ParameterPhy> &pars);
/*!\brief Change the values from a vector of double*/
     virtual void reset(const std::vector<double> &pars);
/*!\brief Rate constant calculator
 *
 * It checks first if the parameters (\f$A\f$ and \f$\beta\f$)
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

/*!\brief \f$\beta\f$ getter*/
     ParameterPhy getBeta()   const {return beta;}
/*!\brief \f$\To\f$ getter*/
     ParameterPhy getT0()     const {return T0;}

/*!\brief \f$\beta\f$ setter*/
     void setBeta(const ParameterPhy &pow) {beta.replace(pow);betaCal.replace(pow);}
/*!\brief \f$\beta\f$ setter*/
     void resetBeta(double pow)            {resetPar(beta,pow);resetPar(betaCal,pow);}
/*!\brief \f$\To\f$ setter*/
     void setT0(const ParameterPhy &Tref)  {T0.replace(Tref);}

/*!\brief Parameters getter*/
     const std::vector<ParameterPhy> getParameters() const;
/*!\brief Number of parameters getter*/
     virtual unsigned int nPars()                    const {return 2;}
/*!\brief Number of optionnal parameters getter*/
     virtual unsigned int nOptionPars()              const {return 1;}

/*!\brief Default configuration: unit*/
     const std::string getDefaultBetaUnit()          const {return std::string();}
/*!\brief Default configuration: min*/
     double getDefaultBetaMin()                      const {return 0.;}
/*!\brief Default configuration: max*/
     double getDefaultBetaMax()                      const {return 2.;}
/*!\brief Default configuration: all in ParameterPhy*/
     const ParameterPhy getDefaultBeta()             const;

/*!\brief Default configuration: unit*/
     virtual const std::string getDefaultUnit(int i)       const;
/*!\brief Default configuration: min*/
     virtual double getDefaultMin(int i)                   const;
/*!\brief Default configuration: max*/
     virtual double getDefaultMax(int i)                   const;
/*!\brief Default configuration: all in ParameterPhy*/
     virtual const ParameterPhy getDefaultParameter(int i) const;

/*!\brief Assignement operator*/
     HercourtEssen &operator=(const HercourtEssen &rhs);

   protected:
/*\brief \f$\beta\f$ parameter*/
     ParameterPhy beta;
/*\brief \f$\beta\f$ parameter used for calculations*/
     ParameterPhy betaCal;
/*\brief \f$\To\f$ parameter*/
     ParameterPhy T0;

/*\brief Check all parameters*/
     virtual void checkParams() const;
/*\brief Sort a vector of parameters
 *
 * We find in this order:
 *   - look for a non-united parameter, this is beta;
 *   - if found, the parameter with exactly K as unit is
 *     the reference temperature, if not found, is it set as
 *     the default (300 K);
 *   - the remaining parameter is the pre-exponential factor.
 */
     virtual std::vector<const ParameterPhy*> whatParametersAreWe(const std::vector<ParameterPhy> &pars);

/*!\brief showKinModel()*/
     void showKinModel(std::ostream &out = std::cout) const;

    private:
/*\brief Templated equation
 *
 * This is the only place where the equation is written,
 * used in the ParameterPhy and double version.
 */
     template<class M>
     M kEquation(const M& temp, const M& pre, const M& p, const M& tref) const;
};
}
#endif
