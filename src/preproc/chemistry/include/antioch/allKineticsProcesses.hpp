#ifndef _ALL_PROCESSES_ARE_HERE_
#define _ALL_PROCESSES_ARE_HERE_

#include "antioch/ElementaryProcess.hpp"
#include "antioch/DuplicateProcess.hpp"
#include "antioch/ThreeBody.hpp"
#include "antioch/LindemannFalloff.hpp"
#include "antioch/TroeFalloff.hpp"

namespace Antioch{
namespace ChemicalProcesses{

enum CHEMICAL_PROCESSES{
  ELEMENTARY_PROCESS,
  DUPLICATE_PROCESS,
  THREE_BODY,
  LINDEMANN_FALLOFF,
  TROE_FALLOFF
};

const std::string ELEMENTARY_PROCESS_STR("Elementary process");
const std::string DUPLICATE_PROCESS_STR("Duplicate process");
const std::string THREE_BODY_STR("Three body");
const std::string LINDEMANN_FALLOFF_STR("Lindemann falloff");
const std::string TROE_FALLOFF_STR("Troe falloff");

const std::string CHEMICAL_PROCESSES_STR[] = {
        ELEMENTARY_PROCESS_STR,
        DUPLICATE_PROCESS_STR,
        THREE_BODY_STR,
        LINDEMANN_FALLOFF_STR,
        TROE_FALLOFF_STR
};

};
namespace BUCKchemicalProcess{
Process * getChemicalProcess(const std::string &key, const Process * copy = NULL);
template <class T>
T *makeChemicalPtr(T * out, const Process * copy = NULL)
{
  if(copy == NULL)
  {
    out = new T;
  }else
  {
    out = new T(*(static_cast<const T*>(copy)));
  }
  return out;
}
};
}

#endif
