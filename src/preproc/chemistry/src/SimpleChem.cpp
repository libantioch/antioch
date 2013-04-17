//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
#include "antioch/SimpleChem.hpp"

namespace Antioch{
SimpleChem::SimpleChem(const SimpleChem &rhs)
{ 
  *this = rhs;
}

void SimpleChem::setMolecule(std::string mol, int imol)
{
  const std::string method("void SimpleChem::setMolecule(std::string , int)");
  if(imol == 0 && molecules.empty())
  {
    addMolecule(mol);
  }else if(imol > (int)molecules.size() - 1)
  {
     antiochError(method,"Requesting out-of-range molecule.");
  }else
  {
    molecules[imol] = mol;
  }
}


int SimpleChem::IsInMolecules(const std::string &mol) const
{
  int out(-1);
  for(unsigned int k = 0; k < molecules.size(); k++)
  {
    if(molecules[k] == mol)
    {
       out = (int)k;
       break;
    }
  }

  return out;
}

bool SimpleChem::areSameMolecules(const std::vector<std::string> &mols) const 
{
  std::vector<unsigned int> skip;
  bool out(true);
  for(unsigned int i = 0; i < mols.size(); i++)
  {
    bool is(false);
    for(unsigned int j = 0; j < molecules.size(); j++)
    {
      bool in(false);
      for(unsigned int k = 0; k < skip.size(); k++)
      {
        in = (j == skip[k]);
      }
      if(in)continue;
      is = (is || (molecules[j] == mols[i]));
      if(is)
      {
        skip.push_back(j);
        break;
      }
    }
    out = (out && is);
  }
  return out;
}

void SimpleChem::eraseMolecule(const std::string &mol)
{
  const std::string method("void SimpleChem::eraseMolecule(const std::string&)");
  int ind(IsInMolecules(mol));
  if(ind == -1)
  {
    antiochWarning(method,"The molecule \"" + mol + "\" is not in the SimpleChem object. Cannot erase it.");
  }else
  {
    eraseMolecule(ind);
  }
  
}

void SimpleChem::setParameter(const ParameterPhy &par, int ipar)
{
  const std::string method("void SimpleChem::setParameter(const ParameterPhy&, int)");
  if(ipar > getNPar() - 1)
        antiochError(method,"Requested parameter does not exist, its index is over the last parameter.");

  parameter[ipar] = par;
}

void SimpleChem::setParameter(const ParameterPhy &par, std::string namepar)
{
  if(namepar.empty())namepar = par.namePar();
  int ind = getIndexByName(namepar,"PP");
  (ind != -1)?setParameter(par,ind):addParameter(par);
}

void SimpleChem::setParameters(const std::vector<ParameterPhy> &pars)
{
  clearParameters();
  for(unsigned int i = 0; i < pars.size(); i++)
  {
    parameter.push_back(pars[i]);
  }
}

int SimpleChem::getIndexByName(const std::string &name,const std::string &type) const
{
  const std::string method("int SimpleChem::getIndexByName(const std::string &,const std::string &) const");

  if(type == "PP")
  {
    for(unsigned int i = 0; i < parameter.size() ; i++)
    {
      if(parameter[i].namePar() == name)return i; 
    }
  }else if(type == "PF")
  {
    for(unsigned int i = 0; i < function.size() ; i++)
    {
      if(function[i].name == name)return i;
    }
  }else if(type == "PHF")
  {
    for(unsigned int i = 0; i < harmFunct.size() ; i++)
    {
      if(harmFunct[i].getName() == name)return i; 
    }
  }else
  {
     antiochError(method,"What is this type you want to search \"" + type + "\"?");
  }

  return -1;
}

void SimpleChem::addEmptyParameter()
{
  parameter.push_back(ParameterPhy());
}

void SimpleChem::addParameters(const std::vector<ParameterPhy> &Vpar)
{
  for(unsigned int i = 0; i < Vpar.size(); i++)
  {
    parameter.push_back(Vpar[i]);
  }
}

void SimpleChem::addParameterByName(std::string parName)
{
  ParameterPhy dummy(parName);
  parameter.push_back(dummy);
}

const std::vector<std::string> SimpleChem::getParametersName() const
{
  std::vector<std::string> names;
  for(unsigned int i = 0; i < parameter.size(); i++)
  {
    names.push_back(parameter[i].namePar());
  }
  return names;
}

void SimpleChem::setParametersName(const std::vector<std::string> &names)
{
  for(unsigned int i =0; i < parameter.size(); i++)
  {
    setParameterName(names[i],i);
  }
}

const ParameterPhy SimpleChem::getParameter(std::string namePar)const 
{
  const std::string method("const ParameterPhy SimpleChem::getParameter(std::string ) const");
  int ipar = getIndexByName(namePar,"PP");
  if(ipar == -1)
  {
     std::string simpleChemNameDesc("\nDescription: \"");
     simpleChemNameDesc += description + "\"\n";
     simpleChemNameDesc += "Molecule(s): ";
     for(unsigned int i = 0; i < molecules.size(); i++)
     {
        simpleChemNameDesc += "\"" + molecules[i] + "\" ";
     }
     simpleChemNameDesc += "\nParameter(s): ";
     for(unsigned int i = 0; i < parameter.size(); i++)
     {
        simpleChemNameDesc += "\"" + parameter[i].namePar() + "\" ";
     }
     antiochError(method,"Parameter \"" + namePar + 
                "\" not found in SimpleChem object" + simpleChemNameDesc);
     return ParameterPhy();
  }
  return parameter[ipar];
}

void SimpleChem::eraseParameter(const std::string &namePar)
{
  const std::string method("void SimpleChem::eraseParameter(const std::string&)");
  int ipar = getIndexByName(namePar,"PP");
  if(ipar == -1)
  {
     antiochWarning(method,"Parameter \"" + namePar + "\" not found, cannot erase it.");
  }
  parameter.erase(parameter.begin() + ipar - 1);
}

void SimpleChem::setUnitObjectParameter(std::string namePar,Units *target)
{
  const std::string method("void SimpleChem::setUnitObjectParameter(std::string ,Units *)");
  int ipar = getIndexByName(namePar,"PP");
  if(ipar == -1)
  {
     antiochError(method,"Parameter \"" + namePar + "\" not found");
     return;
  }
  parameter[ipar].setUnitObject(target);
}

void SimpleChem::setUnitsObjectParameter(std::string namePar)
{
  const std::string method("void SimpleChem::setUnitsObjectParameter(std::string )");
  int ipar = getIndexByName(namePar,"PP");
  if(ipar == -1)
  {
     antiochError(method,"Parameter \"" + namePar + "\" not found");
     return;
  }
}

void SimpleChem::setParameterUnit(std::string namePar,std::string parUnit)
{
  const std::string method("void SimpleChem::setParameterUnit(std::string ,std::string )");
  int ipar = getIndexByName(namePar,"PP");
  if(ipar == -1)
  {
     antiochError(method,"Parameter \"" + namePar + "\" not found");
     return;
  }
  parameter[ipar].setUnitPar(parUnit);
}

void SimpleChem::changeParameterToSomeUnit(std::string target,std::string namePar)
{
  const std::string method("void SimpleChem::changeParameterToSomeUnit(std::string ,std::string )");
  int ipar = getIndexByName(namePar,"PP");
  if(ipar == -1)
  {
     antiochError(method,"Parameter \"" + namePar + "\" not found");
     return;
  }
  parameter[ipar].changeToSomeUnit(target);
}

CoreUnc *SimpleChem::getUncObjectParameterPtr(std::string namePar)
{
  const std::string method("CoreUnc *SimpleChem::getUncObjectParameterPtr(std::string )");
  int ipar = getIndexByName(namePar,"PP");
  if(ipar == -1)
  {
     antiochError(method,"Parameter \"" + namePar + "\" not found");
     return NULL;
  }
  return parameter[ipar].getUncertaintyObject();
}

void SimpleChem::setUncertaintyObjectParameter(std::string namePar,CoreUnc *target)
{
  const std::string method("void SimpleChem::setUncertaintyObjectParameter(std::string ,CoreUnc *)");
  int ipar = getIndexByName(namePar,"PP");
  if(ipar == -1)
  {
     antiochError(method,"Parameter \"" + namePar + "\" not found");
     return;
  }
  parameter[ipar].setUncertaintyObject(target);
}

void SimpleChem::setParameterUncertaintyType(std::string namePar, std::string parUncType)
{
  const std::string method("void SimpleChem::setParameterUncertaintyType(std::string , std::string )");
  int ipar = getIndexByName(namePar,"PP");
  if(ipar == -1)
  {
     antiochError(method,"Parameter \"" + namePar + "\" not found");
     return;
  }
  parameter[ipar].setTypeDPar(parUncType);
}

void SimpleChem::setParameterValue(double value,std::string namePar, int ival)
{
  const std::string method("void SimpleChem::setParameterValue(double, std::string, int )");
  int ipar = getIndexByName(namePar,"PP");
  if(ipar == -1)
  {
     antiochError(method,"Parameter \"" + namePar + "\" not found");
     return;
  }
  parameter[ipar].setValue(value,ival);
}

void SimpleChem::addParameterValue(double value,std::string namePar)
{
  const std::string method("void SimpleChem::addParameterValue(double ,std::string )");
  int ipar = getIndexByName(namePar,"PP");
  if(ipar == -1)
  {
     antiochError(method,"Parameter \"" + namePar + "\" not found");
     return;
  }
  parameter[ipar].addValue(value);
}

void SimpleChem::addParameterValues(const std::vector<double> &vals, std::string namePar)
{
  const std::string method("void SimpleChem::addParameterValues(const std::vector<double> &, std::string )");
  int ipar = getIndexByName(namePar,"PP");
  if(ipar == -1)
  {
     antiochError(method,"Parameter \"" + namePar + "\" not found");
     return;
  }
  parameter[ipar].addValues(vals);
}

void SimpleChem::clearParameter(std::string namePar)
{
  const std::string method("void SimpleChem::clearParameter(std::string )");
  int ipar = getIndexByName(namePar,"PP");
  if(ipar == -1)
  {
     antiochError(method,"Parameter \"" + namePar + "\" not found");
     return;
  }
  parameter[ipar].clear();
}

void SimpleChem::clearParameterValues(std::string namePar)
{
  const std::string method("void SimpleChem::clearParameterValues(std::string )");
  int ipar = getIndexByName(namePar,"PP");
  if(ipar == -1)
  {
     antiochError(method,"Parameter \"" + namePar + "\" not found");
     return;
  }
  parameter[ipar].clearValues();
}

void SimpleChem::setParameterDValue(double value,std::string namePar,int ival)
{
  const std::string method("void SimpleChem::setParameterDValue(double ,std::string ,int )");
  int ipar = getIndexByName(namePar,"PP");
  if(ipar == -1)
  {
     antiochError(method,"Parameter \"" + namePar + "\" not found");
     return;
  }
  parameter[ipar].setDValue(value,ival);
}

void SimpleChem::addParameterDValue(double value,std::string namePar)
{
  const std::string method("void SimpleChem::addParameterDValue(double ,std::string )");
  int ipar = getIndexByName(namePar,"PP");
  if(ipar == -1)
  {
     antiochError(method,"Parameter \"" + namePar + "\" not found");
     return;
  }
  parameter[ipar].addDValue(value);
}

void SimpleChem::addParameterDValues(const std::vector<double> &vals, std::string namePar)
{
  const std::string method("void SimpleChem::addParameterDValues(const std::vector<double> &, std::string )");
  int ipar = getIndexByName(namePar,"PP");
  if(ipar == -1)
  {
     antiochError(method,"Parameter \"" + namePar + "\" not found");
     return;
  }
  parameter[ipar].addDValues(vals);
}

void SimpleChem::resizePar(int newSize,std::string namePar)
{
  const std::string method("void SimpleChem::resizePar(int ,std::string )");
  int ipar = getIndexByName(namePar,"PP");
  if(ipar == -1)
  {
     antiochError(method,"Parameter \"" + namePar + "\" not found");
     return;
  }
  parameter[ipar].resizePar(newSize);
}

void SimpleChem::addParametersByName(const std::vector<std::string> &names)
{
  for(unsigned int i =0; i < names.size(); i++)
  {
    addParameterByName(names[i]);
  }
}

void SimpleChem::addParameterByUnit(std::string parUnit)
{
  ParameterPhy dummy;
  dummy.setUnitPar(parUnit);
  parameter.push_back(dummy);
}


void SimpleChem::densifyParametersUnit()
{
  for(unsigned int i = 1; i < parameter.size(); i++)
  {
    int k = i - 1;
    while(k >= 0)
    {
      if(parameter[i].densifyUnitObject(parameter[k])){break;}
      k--;
    }
  }
}


void SimpleChem::densifyParametersUnitExternal(const SimpleChem &out)
{
  for(unsigned int i = 0; i < parameter.size(); i++)
  {
    int k = (out.getNPar() - 1) + i;
    while(k >= 0)
    {
      if(k < out.getNPar())
      {
        if(parameter[i].densifyUnitObject(out.getParameter(k))){break;}
      }else
      {
        if(parameter[i].densifyUnitObject(parameter[k - out.getNPar()])){break;}
      }
      k--;
    }
  }
  
}

void SimpleChem::densifyParametersToThis(const ParameterPhy &rhs)
{
  for(int ipar = 0; ipar < getNPar(); ipar++)
  {
    parameter[ipar].densifyUnitObject(rhs);
  }
}

void SimpleChem::clearParameters()
{
  for(int i = 0; i < getNPar(); i++)
  {
     clearParameter(i);
  }
}

void SimpleChem::clearParametersValues()
{
  for(int i = 0; i < getNPar(); i++)
  {
     clearParameterValues(i);
  }
}

const std::vector<std::string> SimpleChem::getParametersUnit() const
{
  std::vector<std::string> units;
  for(unsigned int i = 0; i < parameter.size(); i++)
  {
    units.push_back(parameter[i].unitPar());
  }
  return units;
}

void SimpleChem::setParametersUnit(const std::vector<std::string> &units)
{ 
  for(unsigned int i =0; i < parameter.size(); i++)
  {
    setParameterUnit(i,units[i]);
  }
}

void SimpleChem::addParametersByUnit(const std::vector<std::string> &units)
{
  for(unsigned int i =0; i < units.size(); i++)
  {
    addParameterByUnit(units[i]);
  }
}

void SimpleChem::addParameterByUncertaintyType(std::string parUncType)
{
  ParameterPhy dummy;
  dummy.setTypeDPar(parUncType);
  parameter.push_back(dummy);
}

const std::vector<std::string> SimpleChem::getParametersUncertaintyType() const
{
  std::vector<std::string> UncTypes;
  for(unsigned int i = 0; i < parameter.size(); i++)
  {
    UncTypes.push_back(parameter[i].typeDPar());
  }
  return UncTypes;
}

void SimpleChem::setParametersUncertaintyType(const std::vector<std::string> &parUncType)
{
  for(unsigned int i = 0; i < parameter.size(); i++)
  {
    parameter[i].setTypeDPar(parUncType[i]);
  }
}

void SimpleChem::addParametersByUncertaintyType(const std::vector<std::string> &parUncType)
{
  for(unsigned int i = 0; i < parUncType.size(); i++)
  {
    addParameterByUncertaintyType(parUncType[i]);
  }
}

ParameterPhy & SimpleChem::getParameter(const std::string & namePar)
{
  const std::string method("ParameterPhy &SimpleChem::getParameter(std::string )");
  int ipar = getIndexByName(namePar,"PP");
  if(ipar == -1)
  {
     antiochError(method,"Parameter \"" + namePar + "\" not found");
     return parameter[0];
  }
  return parameter[ipar];
}

const void SimpleChem::showAll(std::ostream &out)const 
{
  out << "\n\n************************************************" << std::endl;
  out << "Full description of a SimpleChem" << std::endl;
  out << "Molecule(s) involved: \t";
  for(unsigned int i = 0; i < molecules.size(); i++){out << molecules[i] << ", ";}
  out << std::endl;
  out << "Description: \t" << description << std::endl;
  out << "Number of parameter(s), function(s), multifunction(s): " << parameter.size() << ", " << function.size() << ", " << harmFunct.size() << std::endl;
  if(parameter.size() != 0)
  {
    out << "Parameter(s):\n";
    for(unsigned int i = 0; i < parameter.size(); i++)
    {
      parameter[i].showAll(out);
    }
    out << std::endl;
  }
  if(function.size() != 0)
  {
    out << "Function(s) ({name|unit|uncertainty type}_x,y):\n";
    for(unsigned int i = 0; i < function.size(); i++)
    {
      function[i].showAll(out);
    }
    out << std::endl;
  }
  if(harmFunct.size() != 0)
  {
    out << "Multifunction(s) ({name|unit|uncertainty type}_nOrd):\n";
    for(unsigned int i = 0; i < harmFunct.size(); i++)
    {
      harmFunct[i].showAll(out);
    }
    out << std::endl;
  }
  out << std::endl;
  out << "************************************************\n\n" << std::endl;
}

void SimpleChem::complete(const SimpleChem &rhs)
{
//first, the molecule(s)
  bool isIn = false;
  for(int i = 0; i < rhs.getNMolecules(); i++)
  {
    isIn = false;
    for(unsigned int j = 0; j < molecules.size(); j++)
    {
       isIn = (isIn || (molecules[j] == rhs.getMolecule(i)));
    }
    if(!isIn)molecules.push_back(rhs.getMolecule(i));
  }
//then the description
  if(description != rhs.getDescription())description += rhs.getDescription();

//then the parameter(s)
  for(int i = 0; i < rhs.getNPar(); i++)
  {
    isIn = false;
    for(unsigned int j = 0; j < parameter.size(); j++)
    {
      isIn = (isIn || (parameter[j] == rhs.getParameter(i))); //exists ?
      if(parameter[j] == rhs.getParameter(i))
      {
        parameter[j].complete(rhs.getParameter(i));//if yes, adding info
        break;
      }
    }
    if(!isIn)parameter.push_back(rhs.getParameter(i));//if not, adding parameter
  }

//then the function(s)
//compare the names, add if not in, no order considerations (useless)
  for(int i = 0; i < rhs.getNFunct(); i++)
  {
    bool yep(true);
    for(int j = 0; getNFunct(); j++)
    {
      if(rhs.getFunction(i).name == getFunction(j).name)
      {
        yep = false;
        break;
      }
    }
    if(yep)addFunction(rhs.getFunction(i));
  }

//then the harmFunct(s)
  for(int i = 0; i < rhs.getNHarmFunct(); i++)
  {
     bool yep(true);
     for(int j = 0; j < getNHarmFunct(); j++)
     {
       if(harmFunct[j] == rhs.getHarmFunct(i))
       {
         yep = false;
         break;
       }
     }
    if(yep)addHarmFunct(rhs.getHarmFunct(i));
  }

}

void SimpleChem::setFunction(const PhyFunct &func,int i)
{
  function[i].x = func.x;
  function[i].y = func.y;
  function[i].name = func.name;
}

void SimpleChem::addFunctions(const std::vector<PhyFunct> &funcs)
{
  for(unsigned int i = 0; i < funcs.size(); i++)
  {
     function.push_back(funcs[i]);
  }
}

void SimpleChem::addHarmFuncts(const std::vector<PhyHarmFunct> &relas)
{
  for(unsigned int i = 0; i < relas.size(); i++)
  {
     harmFunct.push_back(relas[i]);
  }
}

PhyFunct SimpleChem::getFunction(const std::string &name) const
{
  const std::string method("PhyFunct Simplechem::getFunction(const std::string &) const");
  int ind = getIndexByName(name,"PF");
  if(ind == -1)antiochError(method,"\"" + name + "\" function is not found in SimpleChem " + molecules[0] + ".");

  return getFunction(ind);
}

PhyHarmFunct SimpleChem::getHarmFunct(const std::string &name) const
{
  const std::string method("PhyHarmFunct SimpleChem::getHarmFunct(const std::string &) const");
  int ind = getIndexByName(name,"PHF");
  if(ind == -1)antiochError(method,name + " multifunction is not found.");

  return getHarmFunct(ind);
}


void SimpleChem::clear()
{
  clearDescription();
  clearMolecules();
  eraseParameters();
}

void SimpleChem::equalizeSC(const SimpleChem &rhs)
{
  molecules = rhs.getMolecules();
  parameter = rhs.getParameters();
  function = rhs.getFunctions();
  harmFunct = rhs.getHarmFuncts();
  description = rhs.getDescription();
}

SimpleChem & SimpleChem::operator = (const SimpleChem &rhs)
{
  if(this == &rhs)return *this;
  equalizeSC(rhs);
  return *this;
}
}
