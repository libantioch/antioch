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
#ifndef _DUPLICATE_CHEMICAL_PROCESS_
#define _DUPLICATE_CHEMICAL_PROCESS_
#include "antioch/Process.hpp"

namespace Antioch{
/*!\file DuplicateProcess.hpp
 * \brief Contains the corresponding class
 *
 * \class DuplicateProcess
 * \brief The simplest chemical reaction
 *
 * It is basically a Process object with
 * a default of Kooij kinetics.
 */

class DuplicateProcess:public Process
{
  public:
/*!\brief Default constructor*/
    DuplicateProcess(const std::string &mod = "Kooij",int nProcess = 2);
/*!\brief Copy constructor*/
    DuplicateProcess(const DuplicateProcess &rhs);
/*!\brief Constructing constructor*/
    DuplicateProcess(const std::string &kinMod, const std::vector<std::vector<ParameterPhy> > &pars);
/*!\brief Default destructor*/
    ~DuplicateProcess();

/*!\brief Initialization
 *
 * 1 - set kinetics model
 * 2 - For all the given processes (size of pars), set the rate constant
 */
        void init(const std::string &kinMod, const std::vector<std::vector<ParameterPhy> > &pars);
/*!\brief Initialization
 *
 * From a vector to a vector of vectors
 */
        void init(const std::string &kinMod, const std::vector<ParameterPhy> &pars) {init(kinMod,sortParameters(pars));}

/*!\brief Rate constant calculator
 *
 * All is defined, all rates constant:
 *   - T is set with the pointer
 */
        virtual std::vector<ParameterPhy> allRateConstant() const;

/*!\brief Rate constant calculator
 *
 *  All is defined internally:
 *    - T
 */
        ParameterPhy rateConstant(int nk = 0) const;
/*!\brief Rate constant calculator
 *
 * Temperature is provided.
 */
        double getRateConstantT(double T, int i = 0) const;

/*!\brief Sort the parameters*/
        std::vector<std::vector<ParameterPhy> > sortParameters(const std::vector<ParameterPhy> &pars);
/*!\brief Sets the temperature*/
        void setTemperature(ParameterPhy *T, int nT = 0);
/*!\brief Parameters setters
 * 
 * If the rate constant exists, it changes the values, if not, it
 * creates it.
 */
        void setAllParameters(const std::vector<std::vector<ParameterPhy> > &pars);
/*!\brief All parameters getter*/
        const std::vector<ParameterPhy> getParameters() const;
/*!\brief All parameters, sorted, getter*/
        const std::vector<std::vector<ParameterPhy> > getAllParameters() const;

};
}
#endif
