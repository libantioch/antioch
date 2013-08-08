#include "antioch/kineticsCrossRoad.hpp"

namespace Antioch{
KineticsModel *kineticsModeling::modelChoice(const std::string &model,const std::vector<ParameterPhy> &pars)
{
  if(model == kineticsModeling::kineticsModel[0] ||
     model == kineticsModeling::kineticsModel[1])
  {
    HercourtEssen *chosenModel = (pars.empty())?new HercourtEssen():new HercourtEssen(pars);
    return static_cast<KineticsModel*>(chosenModel);

  }else if(model == kineticsModeling::kineticsModel[2])
  {
    Arrhenius *chosenModel = (pars.empty())?new Arrhenius():new Arrhenius(pars);
    return static_cast<KineticsModel*>(chosenModel);

  }else if(model == kineticsModeling::kineticsModel[3])
  {
    Berthelot *chosenModel = (pars.empty())?new Berthelot():new Berthelot(pars);
    return static_cast<KineticsModel*>(chosenModel);

  }else if(model == kineticsModeling::kineticsModel[4])
  {
    Kooij *chosenModel = (pars.empty())?new Kooij():new Kooij(pars);
    return static_cast<KineticsModel*>(chosenModel);

  }else if(model == kineticsModeling::kineticsModel[5] ||
           model == kineticsModeling::kineticsModel[6])
  {
    MultHercourtEssen *chosenModel = (pars.empty())?new MultHercourtEssen():new MultHercourtEssen(pars);
    return static_cast<KineticsModel*>(chosenModel);

  }else if(model == kineticsModeling::kineticsModel[7] ||
           model == kineticsModeling::kineticsModel[8])
  {
    BerthelotHercourtEssen *chosenModel = (pars.empty())?new BerthelotHercourtEssen():new BerthelotHercourtEssen(pars);
    return static_cast<KineticsModel*>(chosenModel);

  }else if(model == kineticsModeling::kineticsModel[9] ||
           model == kineticsModeling::kineticsModel[10])
  {
    VantHoff *chosenModel = (pars.empty())?new VantHoff():new VantHoff(pars);
    return static_cast<KineticsModel*>(chosenModel);

  }else
  {
    const std::string method("KineticsModel *ModelChoice(const std::string &,const std::vector<ParameterPhy> &)");
    std::string errorMess = "Your requested model \"" + model + "\" is either not implemented yet or contains a typo.\n";
    errorMess += "Supported kinetics models are:\n";
    for(int i = 0; i < kineticsModeling::nModels; i++)
    {
        errorMess += "\"" + kineticsModeling::kineticsModel[i] + "\", ";
        if(i == 1 || //just so it's more readable
           i == 2 || 
           i == 3 ||
           i == 4 ||
           i == 6 ||
           i == 8 ||
           i == 10)errorMess += "\n";
    }
    antiochError(method,errorMess);
  }

  return NULL;
}
}
