#include "antioch/allKineticsProcesses.hpp"
namespace Antioch{
Process * BUCKchemicalProcess::getChemicalProcess(const std::string &key, const Process * copy)
{
  const std::string method("Process * BUCKchemicalProcess::getChemicalProcess(const std::string &, const Process * (= NULL))");

  if(key == ChemicalProcesses::CHEMICAL_PROCESSES_STR[(int)ChemicalProcesses::ELEMENTARY_PROCESS])
  {
    ElementaryProcess * EP(NULL);
    return static_cast<Process*> (makeChemicalPtr<ElementaryProcess>(EP,copy));
  }else if(key == ChemicalProcesses::CHEMICAL_PROCESSES_STR[(int)ChemicalProcesses::DUPLICATE_PROCESS])
  {
    DuplicateProcess * DP(NULL);
    return static_cast<Process*> (makeChemicalPtr<DuplicateProcess>(DP,copy));
  }else if(key == ChemicalProcesses::CHEMICAL_PROCESSES_STR[(int)ChemicalProcesses::THREE_BODY])
  {
    ThreeBody * TB(NULL);
    return static_cast<Process*> (makeChemicalPtr<ThreeBody>(TB,copy));
  }else if(key == ChemicalProcesses::CHEMICAL_PROCESSES_STR[(int)ChemicalProcesses::LINDEMANN_FALLOFF])
  {
    LindemannFalloff * LF(NULL);
    return static_cast<Process*> (makeChemicalPtr<LindemannFalloff>(LF,copy));
  }else if(key == ChemicalProcesses::CHEMICAL_PROCESSES_STR[(int)ChemicalProcesses::TROE_FALLOFF])
  {
    TroeFalloff * TF(NULL);
    return static_cast<Process*> (makeChemicalPtr<TroeFalloff>(TF,copy));
  }else
  {
    std::string errStr("The asked kinetics model \"" + key + "\" is not supported (yet). Supported models are\n");
    errStr += ChemicalProcesses::CHEMICAL_PROCESSES_STR[(int)ChemicalProcesses::ELEMENTARY_PROCESS] + ", ";
    errStr += ChemicalProcesses::CHEMICAL_PROCESSES_STR[(int)ChemicalProcesses::DUPLICATE_PROCESS]  + ", ";
    errStr += ChemicalProcesses::CHEMICAL_PROCESSES_STR[(int)ChemicalProcesses::THREE_BODY]  + ", ";
    errStr += ChemicalProcesses::CHEMICAL_PROCESSES_STR[(int)ChemicalProcesses::LINDEMANN_FALLOFF]  + " and ";
    errStr += ChemicalProcesses::CHEMICAL_PROCESSES_STR[(int)ChemicalProcesses::TROE_FALLOFF]       + ".";
    antiochError(method,errStr);
  }

  return NULL;
}
}
