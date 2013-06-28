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
#ifndef _KINETICS_CHEMISTRY_
#define _KINETICS_CHEMISTRY_

#include "antioch/SimpleChem.hpp"
#include "antioch/allKineticsProcesses.hpp"
namespace Antioch{

/*!\file KineticsChem.hpp
 * \brief Contains the class KineticsChem
 *
 *\class KineticsChem
 * \brief A reactive SimpleChem, SimpleChem with a rate constant
 *
 */

class KineticsChem:public SimpleChem
{
  public:
    KineticsChem():
        kinProcess(ChemicalProcesses::CHEMICAL_PROCESSES_STR[(int)ChemicalProcesses::ELEMENTARY_PROCESS]),
        reactionRate(NULL),storage(NULL){}
    KineticsChem(const KineticsChem &rhs):reactionRate(NULL),storage(NULL) {*this = rhs;}
    ~KineticsChem();

    void initialize(const SimpleChem &SC, const std::string &kinProc, const Process * kinObj);
    void initialize(const SimpleChem &SC, const std::string &kinProc, const std::string &kinMod = "Kooij", //empty k
                    const std::vector<ParameterPhy> &pars = std::vector<ParameterPhy>(), // EP or DP
                    const std::vector<ParameterPhy> &parsinf = std::vector<ParameterPhy>()); // falloff

    void setNonReactivePart(const SimpleChem &SC)    {equalizeSC(SC);}
    void setKineticsOnly   (const KineticsChem &KC);
    void setKinetics       (const Process *kinObj);
    void setKinetics       (const std::vector<ParameterPhy> &pars = std::vector<ParameterPhy>(), //EP or DP
                            const std::vector<ParameterPhy> &parsinf = std::vector<ParameterPhy>()); // falloff

    void setKineticsProcess(const std::string &kinP) {kinProcess = kinP;}
    void setKineticsModel  (const std::string &kinM) {kinModel   = kinM;}
    void setReactionRate   (const std::string &kinMod = "Kooij", 
                            const Process * kinObj = NULL);

    void storeParameters(const std::vector<ParameterPhy> &pars);
    void storeParameter (const ParameterPhy &par);
    void setReactionFromStorage();
    void setKineticsTemperature  (ParameterPhy *temp);
    void setKineticsConcentration(ParameterPhy *conc);

    const SimpleChem getSimpleChemPart()                    const {return *(static_cast<const SimpleChem*>(this));}
    const std::vector<ParameterPhy> getKineticsParameters() const;

    const std::string getKineticsProcess() const {return kinProcess;}
    const std::string getKineticsModel()   const;

    int nKinPars()                         const;

    std::vector<ParameterPhy> allRateConstant() const;
    ParameterPhy rateConstant(int nk = 0) const;
    double getRateConstantT(double t, int i = 0) const;

//falloff only
    ParameterPhy rateConstantT(double t, int i = 0);
    ParameterPhy rateConstantM(double m, int i = 0);
    double getRateConstantTM(double t, double m, int i = 0) const;

/*!\brief Show all*/
    void showKinetics(std::ostream &out = std::cout) const;


/*!\brief Assignement operator.*/
    KineticsChem & operator=(const KineticsChem &rhs);
/*!\brief Alternative to KineticsChem &operator=(const KineticsChem&)*/
    void equalize(const KineticsChem &rhs) {*this = rhs;}

/*!\brief Process rate getter*/
    const Process * getReaction() const {return reactionRate;}
/*!\brief Parse from a full vector of parameters*/
    void  parseFalloffParameters(std::vector<ParameterPhy> &k0, std::vector<ParameterPhy> &kinf, 
                                 std::vector<ParameterPhy> &Fpars); //if ever other than Troe needs it
    void  parseFalloffParameters(std::vector<ParameterPhy> &k0, std::vector<ParameterPhy> &kinf);

/*!\brief Check if a chemical process is set*/
    bool haveProcess()      const {return !kinProcess.empty();}
/*!\brief Check if a kinetics model is set*/
    bool haveModel()        const {return !kinModel.empty();}
/*!\brief Check if a reaction rate is associated*/
    bool haveReactionRate() const {return (reactionRate != NULL);}

/*!\brief Check if storage is filled*/
    bool haveStorage()      const {return (storage != NULL);}
/*!\brief How many is storage*/
    int howManyStored()     const {return (storage == NULL)?0:(int)storage->size();}

  private:
    std::string kinProcess;
    std::string kinModel;
    Process *reactionRate;

    std::vector<ParameterPhy> *storage;//should be empty

    bool isFalloff() const;
};

}
#endif
