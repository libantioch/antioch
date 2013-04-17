//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
#ifndef _BERTHELOT_KINETICS_MODELS_
#define _BERTHELOT_KINETICS_MODELS_

#include "antioch/KineticsModel.hpp"
#include "antioch/phymath.hpp"

namespace Antioch{
 /*!\file Berthelot.hpp
 * \brief contains corresponding class
 *
 * \class Berthelot
 * \ingroup kinmod
 * \brief The Berthelot kinetics model
 *
 * The Berthelot equation is:
 * \f[
 *     \BerthEquation
 * \f]
 *
 * The uncertainty following the GUM recommandation is:
 * \f[
 * \renewcommand{\arraystretch}{1.5}
 * \begin{array}{rll}
 *    u_{\Berth}^2 =& \left(\BerthEq\right)^2\times \tabularnewline
 *    ( \tabularnewline
 *                  &   \left(\frac{1}{A}\right)^2 & u_{A}^2   \tabularnewline
 *                  &+  T^2                        & u_{D}^2 \tabularnewline
 *                  &+  D^2 & u_{T}^2   \tabularnewline
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
 */
class Berthelot:virtual public KineticsModel
{
   public:
/*!\brief Default constructor*/
     Berthelot(const ParameterPhy &d = ParameterPhy("Exponential parameter","K-1")):D(d),DCal(d){}
/*!\brief Constructor*/
     Berthelot(const ParameterPhy &Pre, const ParameterPhy &d):
        KineticsModel(Pre),D(d),DCal(d){}
/*!\brief Copy constructor*/
     Berthelot(const Berthelot &rhs) {*this = rhs;}
/*!\brief Automatic choosing constructor*/
     Berthelot(const std::vector<ParameterPhy> &pars);
/*!\brief Destructor*/
    virtual ~Berthelot(){}

/*!\brief Initialize from a vector of parameters*/
     virtual void init(const std::vector<ParameterPhy> &pars);
/*!\brief Change the values from a vector of double*/
     virtual void reset(const std::vector<double> &pars);
/*!\brief Rate constant calculator
 *
 * It checks first if the parameters (\f$A\f$ and \f$D\f$)
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
     virtual double getRateConstantT(double Tem, int i = 0) const;

/*!\brief \f$D\f$ getter*/
     ParameterPhy getD()        const {return D;}

/*!\brief \f$D\f$ setter*/
     void setD(const ParameterPhy &d) {D.replace(d);DCal.replace(d);}
/*!\brief \f$D\f$ setter*/
     void resetD(double d)            {resetPar(D,d);resetPar(DCal,d);} 

/*!\brief Parameters getter*/
     const std::vector<ParameterPhy> getParameters() const;
/*!\brief Number of parameters getter*/
     virtual unsigned int nPars() const {return 2;}

/*!\brief Default configuration: unit*/
     const std::string getDefaultDUnit() const {return "K-1";}
/*!\brief Default configuration: min*/
     double getDefaultDMin()             const {return 0.;};
/*!\brief Default configuration: max*/
     double getDefaultDMax()             const {return 1000.;}
/*!\brief Default configuration: all in ParameterPhy*/
     const ParameterPhy getDefaultD()    const;

/*!\brief Default configuration: unit*/
     virtual const std::string getDefaultUnit(int i)       const;
/*!\brief Default configuration: min*/
     virtual double getDefaultMin(int i)                   const;
/*!\brief Default configuration: max*/
     virtual double getDefaultMax(int i)                   const;
/*!\brief Default configuration: all in ParameterPhy*/
     virtual const ParameterPhy getDefaultParameter(int i) const;


/*!\brief Assignement operator*/
     Berthelot &operator=(const Berthelot &rhs);

   protected:
/*\brief \f$D\f$ parameter*/
     ParameterPhy D;
/*\brief \f$D\f$ parameter used for calculations*/
     ParameterPhy DCal;

/*\brief Check all parameters*/
     virtual void checkParams() const;
/*\brief Sort a vector of parameters
 *
 * The check is made on D, unit is imposed as K-1
 */
     virtual std::vector<const ParameterPhy*> whatParametersAreWe(const std::vector<ParameterPhy> &pars);
   private:
/*\brief Templated equation
 *
 * This is the only place where the equation is written,
 * used in the ParameterPhy and double version.
 */
     template<class M>
     M kEquation(const M& temp, const M& pre, const M& d) const;
};
}
#endif
