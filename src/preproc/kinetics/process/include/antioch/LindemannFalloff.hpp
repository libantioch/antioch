//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#ifndef _LINDEMANN_FALLOFF_
#define _LINDEMANN_FALLOFF_

#include "antioch/FalloffProcess.hpp"

namespace Antioch{
/*!\file LindemannFalloff.hpp
 * \brief Lindemann falloff chemical process
 *
 * This file contains the LindemannFalloff class.
 *
 * \class LindemannFalloff
 * \brief Lindemann falloff class
 *
 * This is the simplest falloff, \f$F = 1\f$.
 * This class adds nothing else to the FalloffProcess
 * base class.
 *
 */

class LindemannFalloff:public FalloffProcess
{
     public:
/*!\brief Default constructor*/
        LindemannFalloff(){setKineticsProcess("Lindemann falloff");}
/*!\brief Copy constructor*/
        LindemannFalloff(const LindemannFalloff &rhs);
/*!\brief Constructor*/
        LindemannFalloff(const std::string &kinMod, const std::vector<ParameterPhy> &pars0, const std::vector<ParameterPhy> &parsinf, ParameterPhy *m = NULL):
             FalloffProcess(kinMod,pars0,parsinf,m){setKineticsProcess("Lindemann falloff");}
             
/*!\brief Destructor*/
        ~LindemannFalloff(){}

    private:
     ParameterPhy FT(double t, int i = 0)      const {return FLind;}
     ParameterPhy FM(double m, int i = 0)      const {return FLind;}
     double FTM(double t, double m, int i = 0) const {return 1.0;}
     /*\brief reset the parameters for F*/
     void resetF(std::vector<double> &pars){}

/*!\brief Virtual pure function for double version of F
 *
 * Here M is defined, T is provided
 */
     ParameterPhy logFT(double t, int i = 0)         {return log(FLind);}
/*!\brief Virtual pure function for double version of F
 *
 * Here T is defined, M is provided
 */
     ParameterPhy logFM(double m, int i = 0)         {return log(FLind);}
/*!\brief Virtual pure function for double version of F
 *
 * Here all (T and M) is provided
 */
     double logFTM(double t, double m, int i = 0) const {return 0.0;}


     static ParameterPhy FLind;
};
}
#endif
