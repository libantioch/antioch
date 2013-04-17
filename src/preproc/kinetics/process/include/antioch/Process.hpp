//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
#ifndef _REACTION_BASE_
#define _REACTION_BASE_

#include "antioch/ParameterPhy.hpp"
#include "antioch/phymath.hpp"
#include "antioch/kineticsCrossRoad.hpp" //include the hpps

namespace Antioch{
/*!\file Process.hpp
 * \brief This file contains the basics for reactions
 *
 * It contains the class Process and the public method
 * to choose the process from the std::string key.
 *
 *\class Process 
 * \ingroup kinproc
 * \brief Basic class usage of all reactions
 *
 * This class contains the simplest operations for a chemical
 * reaction. Thus a kinetics model and a rate constant object.
 * The calculations are done by the rate constant object, given
 * a double temperature of a ParameterPhy adress.
 *
 * \defgroup kinproc Kinetics processes
 * The kinetics processes are the representation of ``how
 * elementary'' the process is. Possibilities are:
 *   - the reaction is indeed an elementary process, thus
 *      we use the ElementaryProcess object;
 *   - the reaction is not well approximated by an
 *      elementary process, thus we use a DuplicateProcess
 *      object.
 *
 *  If the reaction depends on the pressure (\e i.\e e. the
 *  concentration of the medium), we need to consider 
 *  the three-body reactions (class ThreeBody). Three-body reactions are 
 *  basically bimolecular pressure-dependant reactions,
 *  which rate constant is a multiple of the bimolecular
 *  reaction.
 *
 *  Falloff reactions are 
 */
class Process{
      public:
/*!\brief Default constructor*/
        Process():M(NULL){rates.push_back(NULL);}
/*!\brief Copy constructor*/
        Process(const Process &rhs);
/*!\brief constructor for a kinetics model only*/
        Process(const std::string &kinMod):kineticsMod(kinMod),M(NULL){rates.push_back(kineticsModeling::modelChoice(kineticsMod));}
/*!\brief Full constructor*/
        Process(const std::string &kinMod, const std::vector<ParameterPhy> &pars);
/*!\brief Default destructor*/
        virtual ~Process();

/*!\brief Initialization
 *
 * 1 - set kinetics model
 * 2 - set parameters
 */
        virtual void init(const std::string &kinMod, const std::vector<ParameterPhy> &pars);
/*!\brief Initialization for falloff
 *
 * 1 - set kinetics model
 * 2 - set parameters k0
 * 3 - set parameters kinf
 */
        virtual void init(const std::string &kinMod, const std::vector<ParameterPhy> &pars0, const std::vector<ParameterPhy> &parsinf, ParameterPhy *m = NULL){}

/*!\brief Rate constant calculator
 *
 * All is defined:
 *   - T is set with the pointer
 */
        virtual ParameterPhy rateConstant(int nk = 0) const;
/*!\brief Rate constant calculator
 *
 * All is defined, all rates constant:
 *   - T is set with the pointer
 */
        virtual std::vector<ParameterPhy> allRateConstant() const;

/*!\brief Rate constant calculator
 *
 * Here we give the temperature on non falloff
 */
        virtual double getRateConstantT(double T, int i = 0) const;

//******************* falloff ********************//

/*!\brief Rate constant calculator
 *
 * Here we give the temperature on falloff, M is defined
 */
        virtual ParameterPhy rateConstantT(double t, int i = 0) const;
/*!\brief Rate constant calculator
 *
 * Here we give the concentration on falloff, T is defined
 */
        virtual ParameterPhy rateConstantM(double m, int i = 0) const;
/*!\brief Rate constant calculator for falloff
 *
 * All is provided.
 */
        virtual double getRateConstantTM(double t, double m, int i = 0) const;

/*!\brief Sets the temperature
 *
 * The index gives what temperature you want to set:
 *   - 0 (default) is the temperature
 *   - 1 (multiple Hercourt-Essen) is for Ti
 */
        virtual void setTemperature(ParameterPhy *T, int nT = 0);
/*!\brief \f$T\f$ getter*/
        ParameterPhy * getTemperature(int nT = 0) const;
/*!\brief Sets the mixture concentration*/
        virtual void setConcentration(ParameterPhy *m) {M = m;}
/*!\brief Get the mixture concentration*/
        virtual ParameterPhy * getConcentration() const {return M;}
  
/*!\brief Test to see if T is defined*/
        bool isTDefined();
/*!\brief Is the mixture concentration defined*/
        bool isMDefined() {return !(M == NULL);}
/*!\brief Default unit of M*/
        std::string defaultMUnit() const {return "mol/cm3";}

/*!\brief Number of kinetics parameters getter*/
     unsigned int nKinPars(int iProc = 0)       const;
/*!\brief Number of kinetics optionnal parameters getter*/
     unsigned int nKinOptionPars(int iProc = 0) const;
/*!\brief Total number of kinetics parameters getter*/
     unsigned int nTotalKinPars(int iProc = 0)  const;

/* Default prior configuration:
 * - min
 * - max
 * - unit
 *
 *   pdf is uniform => beware preexp, it's logu
 *
 */

/*!\brief Default prior configuration: min value*/
    virtual double getDefaultPriorMinValue(int i)              const;
/*!\brief Default prior configuration: max value*/
    virtual double getDefaultPriorMaxValue(int i)              const;
/*!\brief Default prior configuration: unit*/
    virtual const std::string getDefaultPriorUnit(int i)       const;
/*!\brief Default prior configuration: parameter*/
    virtual const ParameterPhy getDefaultPriorParameter(int i) const;


/*!\brief Parameters setters
 * 
 * If the rate constant exists, it changes the values, if not, it
 * creates it.
 */
        virtual void setParameters(const std::vector<ParameterPhy> &pars, KineticsModel *&rateC);
/*!\brief Parameters resetters
 * 
 * Sends the values to the rate constant.
 */
        void resetParameters(const std::vector<double> &pars);

/*!\brief Kinetics model setter*/
        void setKineticsModel(const std::string &kin)    {kineticsMod = kin;}
/*!\brief Kinetics model getter*/
        const std::string getKineticsModel()       const {return kineticsMod;}
/*!\brief Kinetics process setter*/
        void setKineticsProcess(const std::string &proc) {kineticsProc = proc;}
/*!\brief Kinetics process getter*/
        const std::string getKineticsProcess()     const {return kineticsProc;}

/*!\brief All parameters getter*/
        virtual const std::vector<ParameterPhy> getParameters() const {return getParametersFromRate(rates[0]);}
/*!\brief Parameters getter for one pointer in a vector form*/
        const std::vector<ParameterPhy> getParametersFromRate(KineticsModel * rateC) const;

/*!\brief Sets F parameters, for fallorr*/
        virtual void setFalloffParameters(const std::vector<ParameterPhy> &pars){} //nothing to do here
/*!\brief Number of parameters needed for falloff (F)*/
        virtual unsigned int nFPars()       const {return 0;}
/*!\brief Number of optionnal parameters needed for falloff (F)*/
        virtual unsigned int nFOptionPars() const {return 0;}
/*!\brief Total number of parameters needed for falloff (F)*/
        unsigned int nFTotalPars()          const {return (nFPars() + nFOptionPars());}

/*!\brief Default configuration F: unit*/
        virtual const std::string getDefaultFUnit(int i)            const {return std::string();}
/*!\brief Default configuration F: min*/
        virtual double getDefaultFMin(int i)                        const {return 0.;}
/*!\brief Default configuration F: max*/
        virtual double getDefaultFMax(int i)                        const {return 0.;}
/*!\brief Default configuration F: all in ParameterPhy*/
        virtual const ParameterPhy getDefaultFPriorParameter(int i) const {return ParameterPhy();}
/*!\brief Rate constant ptr getter*/
        const KineticsModel * getRate(int i = 0)                    const {return rates[i];}
/*!\brief Number of kinetics rate*/
        unsigned int getKineticsDegree()                            const {return rates.size();}

/*\brief Numbre of processes getter*/
        unsigned int getNProcesses()                                const {return rates.size();}

/*!\brief showAll()*/
        void showAll(std::ostream &out = std::cout) const;

     protected:
        std::string kineticsMod;
        std::string kineticsProc;
        std::vector<KineticsModel*> rates;

        ParameterPhy *M;

/*!\brief Parameter reset*/
        void resetPar(ParameterPhy &par,double d);
};
}
#endif
