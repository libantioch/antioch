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
#ifndef _KINETICS_MODELS_
#define _KINETICS_MODELS_

#include "antioch/ParameterPhy.hpp"
#include "antioch/phymath.hpp"

namespace Antioch{

/*!\file KineticsModel.hpp
 * \brief contains mother class of kinetics models
 *
 * \class KineticsModel
 * \ingroup kinmod
 * \brief mother of the kinetics model
 *
 * This class is an abstract class, the pure virtual methods are
 *   - virtual double getRateConstantT(double Tem, int i = 0) const = 0;
 *   - virtual ParameterPhy rateConstant(int (=0)) = 0;
 *   - virtual void checkParams() const = 0;
 *   - virtual std::vector<const ParameterPhy*> whatParametersAreWe(const std::vector<ParameterPhy> &pars) = 0;
 *   - virtual void init(const std::vector<ParameterPhy*> &pars) = 0;
 *   - virtual const std::vector<ParameterPhy> getParameters() = 0;
 *   - virtual unsigned int nPars() const = 0;
 *   - virtual void reset(const std::vector<double> &pars) = 0;
 *
 *   \defgroup kinmod Kinetics models
 *   The kinetics models used are separable in two
 *   groups:
 *     - the simple models
 *     - the combined models
 *   The simple models are the simplest equations used for 
 *   kinetics, the combined models are combinations (multiplication)
 *   of the former models.
 *
 *   The simple models are:
 *     - the Hercourt-Essen model (class HercourtEssen)
 *     - the Berthelot model  (class Berthelot)
 *     - the Arrhenius model (class Arrhenius)
 *
 *   The combined models are:
 *     - the Kooij model (class Kooij)
 *     - the Berthelot Hercourt-Essen model (class BerthelotHercourtEssen)
 *     - the Van't Hoff model (class VantHoff)
 *
 *   Note that we also have the two temperatures Hercourt-Essen,
 *   coded by the class MultHercourtEssen.
 *
 *   The equations are similar, they all contain a pre-exponential parameter, and
 *   are a combination of different terms:
 *     - a temperature dependence of the form \f$\left(\frac{T}{\mathrm{T_{ref}}}\right)^\beta\f$,
 *     - an exponential factor of the form \f$\exp\left(-\frac{E_a}{T}\right)\f$
 *     - an exponential factor of the form \f$\exp\left(-DT\right)\f$
 *
 *   The simple models contain only one term, the combined models at least two.
 */
class KineticsModel
{
  public:
/*\brief Constructor*/
     KineticsModel(const ParameterPhy &Pre = ParameterPhy("Pre exponential parameter","cm3/mol/s")):A(Pre),ACal(Pre),T(NULL){}
/*\brief Destructor*/
    virtual ~KineticsModel(){}

//********* pure virtual functions *****************//

/*!\brief Initialize from a vector of parameters*/
     virtual void init(const std::vector<ParameterPhy> &pars)      = 0;
/*!\brief Change the values from a vector of double*/
     virtual void reset(const std::vector<double> &pars)           = 0;
/*!\brief Rate constant calculator*/
     virtual const ParameterPhy rateConstant(int nk = 0)           = 0;
/*!\brief Rate constant calculator*/
     virtual double getRateConstantT(double Tem, int i = 0)  const = 0;
/*!\brief Number of parameters getter*/
     virtual unsigned int nPars()                            const = 0;
/*!\brief Parameters getter*/
     virtual const std::vector<ParameterPhy> getParameters() const = 0;

//********* end pure virtual functions *****************//


//************** rate constant ********************//

/*!\brief All rate constant calculator*/
     const std::vector<ParameterPhy> allRateConstant();


//********* Temperature ***********//

/*!\brief \f$T\f$ setter*/
     virtual void setT(ParameterPhy *Tused, int nT = 0);
/*!\brief \f$T\f$ getter*/
     virtual ParameterPhy * getT(int nT = 0)                 const;
/*!\brief Test to see if T is defined*/
     bool isTDefined()                                              {return !(T == NULL);}


//********* Pre exponential parameter ***********//

/*!\brief \f$A\f$ getter*/
     ParameterPhy getPreExp()                 const {return A;}
/*!\brief \f$A\f$ setter*/
     void setPreExp(const ParameterPhy &Pre)        {A.replace(Pre);ACal.replace(Pre);}
/*!\brief \f$A\f$ setter*/
     void resetPreExp(double Pre)                   {resetPar(A,Pre);resetPar(ACal,Pre);}

/*!\brief Default configuration: unit*/
     const std::string getDefaultPreExpUnit() const {return "cm3/mol/s";}
/*!\brief Default configuration: min*/
     double getDefaultPreExpMin()             const {return 1.;}
/*!\brief Default configuration: max*/
     double getDefaultPreExpMax()             const {return 1e20;}
/*!\brief Default configuration: all in ParameterPhy*/
     const ParameterPhy getDefaultPreExp()    const;


/*!\brief Default configuration: unit*/
     virtual const std::string getDefaultUnit(int i)       const {return "cm3/mol/s";}
/*!\brief Default configuration: min*/
     virtual double getDefaultMin(int i)                   const {return 1.;}
/*!\brief Default configuration: max*/
     virtual double getDefaultMax(int i)                   const {return 1e20;}
/*!\brief Default configuration: all in ParameterPhy*/
     virtual const ParameterPhy getDefaultParameter(int i) const {return getDefaultPreExp();}


//************** Optionnal parameters ****************//

/*!\brief Number of optionnal parameters getter*/
     virtual unsigned int nOptionPars() const {return 0;}
/*!\brief Total number of parameters getter*/
     unsigned int nTotalPars()          const {return (nPars() + nOptionPars());}

//*************** show all ********************//

/*!\brief showAll() */
     void showAll(std::ostream &out = std::cout) const;

  protected:

//********* pure virtual functions *****************//

/*\brief Check all parameters*/
     virtual void checkParams() const = 0;
/*\brief Sort a vector of parameters*/
     virtual std::vector<const ParameterPhy*> whatParametersAreWe(const std::vector<ParameterPhy> &pars) = 0;

//********* end pure virtual functions *****************//


//*************** parameters ************//

/*\brief \f$A\f$ parameter*/
     ParameterPhy A;
/*\brief \f$A\f$ parameter used for calculations (1 value)*/
     ParameterPhy ACal;
/*\brief Pointer to \f$T\f$ parameter*/
     ParameterPhy *T;


//************** show kin & T *************//

/*!\brief Show kinetics parameters*/
     virtual void showKinModel(std::ostream &out = std::cout) const {out << "To be implemented" << std::endl;}
/*!\brief Show temperature parameter*/
     void showT(std::ostream &out = std::cout) const;



//**************** check parameters ********************//

/*\brief Check if the A is empty or not.*/
     void checkA(const std::string &object) const;
/*\brief Check if the T pointer is NULL or not.*/
     void checkT(const std::string &object) const;
/*\brief Check if parameter is empty or not.*/
     void checkParam(const ParameterPhy &test, const std::string &namepar, const std::string &object) const;

//**************** reset parameters ********************//

/*!\brief Parameter reset*/
     void resetPar(ParameterPhy &par, double r);


//**************** misc ********************//

/*\brief General two parameters one defined unit sorting
 *
 * Within a vector of two parameters, find the one having exactly
 * the given unit, push it first, then the other.
 */
     std::vector<const ParameterPhy*> findSameUnitInTwoParameters(const std::vector<ParameterPhy> &pars,const std::string &unit);
/*\brief General two parameters one defined unit sorting
 *
 * Within a vector of two parameters, find the one having an
 * homogeneous unit to
 * the given unit, push it first, then the other.
 */
     std::vector<const ParameterPhy*> findHomoUnitInTwoParameters(const std::vector<ParameterPhy> &pars,const std::string &unit);
};
}
#endif
