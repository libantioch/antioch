//-----------------------------------------------------------------------bl-
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
