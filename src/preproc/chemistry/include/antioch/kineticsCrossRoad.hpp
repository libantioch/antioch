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
#ifndef _KINETICS_CROSS_ROADS_
#define _KINETICS_CROSS_ROADS_

#include "antioch/HercourtEssen.hpp"
#include "antioch/Arrhenius.hpp"
#include "antioch/Berthelot.hpp"
#include "antioch/MultHercourtEssen.hpp"
#include "antioch/BerthelotHercourtEssen.hpp"
#include "antioch/Kooij.hpp"
#include "antioch/VantHoff.hpp"

namespace Antioch{
/*!\file kineticsCrossRoad.hpp
 * \brief All kinetics models plus public function to send a chosen kinetics
 * object
 *
 * \namespace kineticsModeling
 * \brief A namespace to store the key to kinetics model
 */

namespace kineticsModeling{
/*!\brief Array to contain the string for each model*/
const std::string kineticsModel[]={
"Hercourt-Essen","HE",            //Hercourt-Essen keys
"Arrhenius",                      //Arrhenius keys
"Berthelot",                      //Berthelot keys
"Kooij",                          //Kooij keys
"mHercourt-Essen","mHE",          //multiple Hercourt-Essen keys
"Berthelot Hercourt-Essen","BHE", //Berthelot Hercourt-Essen keys
"Van't Hoff","VH"                 //Van't Hoff keys
};
const int nModels(11);

/*\!brief Model chooser
 *
 * This is the bad if..else if..else if..
 * implementation. Follow the kineticsModeling::kineticsModel
 * array of string definition.
 */
KineticsModel *modelChoice(const std::string &model,const std::vector<ParameterPhy> &pars = std::vector<ParameterPhy>());
};
}
#endif
