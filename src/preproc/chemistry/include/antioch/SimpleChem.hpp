//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
#ifndef _KINETICS_SIMPLE_STRUCTURE_
#define _KINETICS_SIMPLE_STRUCTURE_

#include "antioch/ParameterPhy.hpp"
#include "antioch/PhyFunct.hpp"
#include "antioch/PhyHarmFunct.hpp"

namespace Antioch{
/*! \class SimpleChem
 * \brief A class to code chemical entities.
 * \todo Adding PhyFunct and PhyHarmFunct capabilities:
 *   - showAll(std::ostream (=std::cout)) to be completed
 *   - clear() to be completed
 *   - erase functions for objects
 *
 * A simple chemical entity is some molecules associated
 * to some physical properties to describe the
 * any wanted phenomenon related to these molecules.
 *
 * The molecules are simply coded as a std::vector of std::strings,
 * the physical properties as a std::vector of ParameterPhy, and
 * a std::string is used to describe the phenomenon coded. 
 */
class SimpleChem{
        public:
/*! \brief Default constructor.*/
         SimpleChem(){}
/*! \brief Copy constructor.*/
         SimpleChem(const SimpleChem &rhs);
/*! \brief Building constructor.*/
         SimpleChem(const std::vector<std::string> &entities, 
                    const std::vector<ParameterPhy> &properties,
                    std::string des):molecules(entities),parameter(properties),description(des){}
/*! \brief Default destructor.*/
         ~SimpleChem(){}

/*!\brief Description getter*/
         const std::string getDescription()    const {return description;}
/*!\brief Description setter*/
         void setDescription(const std::string &des) {description = des;}
/*!\brief Description clearer*/
         void clearDescription()                     {description.clear();}

/*!\brief Molecules getter*/
         const std::vector<std::string> getMolecules()     const {return molecules;}
/*!\brief Molecules setter*/
         void setMolecules(const std::vector<std::string> &mols) {molecules = mols;}
/*!\brief Number of molecules getter*/
         int getNMolecules()                               const {return molecules.size();}
/*!\brief Molecule comparator, check wether the asked molecule is in the SimpleChem object.
 *
 * It returns the index of the molecule if found, -1 if not found.
 */
         int IsInMolecules(const std::string &mol) const;
/*!\brief Molecules comparator*/
         bool areSameMolecules(const std::vector<std::string> &mols) const;
/*!\brief Molecules clearer*/
         void clearMolecules()                                  {molecules.clear();}

/*!\brief Molecule getter*/
         const std::string getMolecule(int imol = 0)      const {return molecules[imol];}
/*!\brief Molecule getter overload with unsigned int*/
         const std::string getMolecule(unsigned int imol) const {return molecules[imol];}
/*!\brief Molecule adder*/
         void addMolecule(std::string mol)                      {molecules.push_back(mol);}
/*!\brief Molecule setter*/
         void setMolecule(std::string mol, int imol = 0);
/*!\brief Molecule setter overload with unsigned int*/
         void setMolecule(std::string mol, unsigned int imol)   {molecules[imol] = mol;}
/*!\brief Molecule eraser*/
         void eraseMolecule(int imol)                           {molecules.erase(molecules.begin() + imol - 1);}
/*!\brief Molecule eraser, unsigned int overload*/
         void eraseMolecule(unsigned int imol)                  {eraseMolecule((int)imol);}
/*!\brief Molecule eraser, name overload*/
         void eraseMolecule(const std::string& mol);

/*!\brief Parameters getter*/
         const std::vector<ParameterPhy> getParameters()      const {return parameter;}
/*!\brief Parameter getter by index*/
         const ParameterPhy getParameter(int ipar)            const {return parameter[ipar];}
/*!\brief Parameter getter by index, unsigned int overload*/
         const ParameterPhy getParameter(unsigned int ipar)   const {return parameter[ipar];}
/*!\brief Parameter getter by name*/
         const ParameterPhy getParameter(std::string namePar) const;
/*!\brief Parameter setter*/
         void setParameter(const ParameterPhy &par, int ipar);
/*!\brief Parameter setter, unsigned int overload*/
         void setParameter(const ParameterPhy &par, unsigned int ipar) {setParameter(par,(int)ipar);}
/*!\brief Parameter setter, name lookout
 *
 * The search of the parameter to be equalized is done
 * by name equality, meaning if the names are the same,
 * then it's ok, whatever the unit, uncertainty and values.
 *
 * There are two ways of setting a parameter. The question to
 * solve first is to know if this is a strict equality (the name
 * are also the same, see 3PH, AtomicParameter &operator=(const AtomicParameter &)
 * for details), or a numerical equality. The difference is
 * basically the difference between updating a parameter and
 * adding/changing a parameter.
 *
 * In case of an update, only the updated parameter is needed,
 * as the name of the parameter to be changed and the input
 * parameter are the same. Thus the std::string namepar is
 * not required.
 *
 * In case of a complete change of the parameter, the names
 * don't match, and the name of the parameter to be changed
 * is required.
 *
 * Once the name ambiguity is solved, the methods searchs for
 * the parameter. If it is found, then it is equalized, if it
 * is not found, the parameter is added.
 */
         void setParameter(const ParameterPhy &par, std::string namepar = std::string());
/*!\brief Parameters setter*/
         void setParameters(const std::vector<ParameterPhy> &pars);
/*!\brief Number of parameters getter*/
         int getNPar()                                const {return parameter.size();}
/*!\brief Parameter adder*/
         void addParameter(const ParameterPhy &par)         {parameter.push_back(par);}
/*!\brief Parameter adder, only the name is given*/
         void addParameterByName(std::string parName);
/*!\brief Parameter adder, only the unit is given*/
         void addParameterByUnit(std::string parUnit);
/*!\brief Parameters adder*/
         void addParameters(const std::vector<ParameterPhy> &Vpar);
/*!\brief An empty parameter adder*/
         void addEmptyParameter();

/*!\brief Name of parameter getter by index*/
         const std::string getParameterName(int ipar)                 const {return parameter[ipar].namePar();}
/*!\brief Name of parameter getter by index, unsigned int overload*/
         const std::string getParameterName(unsigned int ipar)        const {return parameter[ipar].namePar();}
/*!\brief Name of parameter setter by index*/
         void setParameterName(std::string parName,int ipar)                {parameter[ipar].setNamePar(parName);}
/*!\brief Name of parameter setter by index, unsigned int overload*/
         void setParameterName(std::string parName,unsigned int ipar)       {parameter[ipar].setNamePar(parName);}

/*!\brief Parameters name getter*/
         const std::vector<std::string> getParametersName() const;
/*!\brief Parameters name setter*/
         void setParametersName(const std::vector<std::string> &names);
/*!\brief Parameters name adder*/
         void addParametersByName(const std::vector<std::string> &names);

/*!\brief Parameter unit pointer getter by index*/
         Units *getUnitObjectParameterPtr(int ipar)                const {return parameter[ipar].getUnitObject();}
/*!\brief Parameter unit pointer getter by index, unsigned int overload*/
         Units *getUnitObjectParameterPtr(unsigned int ipar)       const {return parameter[ipar].getUnitObject();}
/*!\brief Parameter unit pointer getter by name*/
         Units *getUnitObjectParameterPtr(std::string namePar)     const {return getParameter(namePar).getUnitObject();}

/*!\brief Parameter unit pointer setter by index*/
         void setUnitObjectParameter(int ipar,Units *target)             {parameter[ipar].setUnitObject(target);}
/*!\brief Parameter unit pointer setter by index, unsigned int overload*/
         void setUnitObjectParameter(unsigned int ipar,Units *target)    {parameter[ipar].setUnitObject(target);}
/*!\brief Parameter unit pointer setter by name*/
         void setUnitObjectParameter(std::string namePar,Units *target);
/*!\brief Parameters setter to Units object unit mode by name*/
         void setUnitsObjectParameter(std::string namePar);

/*!\brief Parameter unit getter by index*/
         const std::string getParameterUnit(int ipar)            const {return parameter[ipar].unitPar();}
/*!\brief Parameter unit getter by index, unsigned int overload*/
         const std::string getParameterUnit(unsigned int ipar)   const {return parameter[ipar].unitPar();}
/*!\brief Parameter unit getter by name*/
         const std::string getParameterUnit(std::string namePar) const {return getParameter(namePar).unitPar();}
/*!\brief Parameter unit setter by index*/
         void setParameterUnit(int ipar, std::string parUnit)          {parameter[ipar].setUnitPar(parUnit);}
/*!\brief Parameter unit setter by index, unsigned int overload*/
         void setParameterUnit(unsigned int ipar, std::string parUnit) {parameter[ipar].setUnitPar(parUnit);}
/*!\brief Parameter unit setter by name*/
         void setParameterUnit(std::string namePar,std::string parUnit);
/*!\brief Parameter factor to SI unit by index
 *
 * This method must be used with an advanced unit object
 */
         double getFactorToSI(int ipar)                    const {return parameter[ipar].getFactor();}
/*!\brief Parameter factor to SI unit by index, unsigned int overload
 *
 * This method must be used with an advanced unit object
 */
         double getFactorToSI(unsigned int ipar)           const {return parameter[ipar].getFactor();}
/*!\brief Parameter factor to SI unit by name
 *
 * This method must be used with an advanced unit object
 */
         double getFactorToSI(std::string namePar)              const {return getParameter(namePar).getFactor();}
/*!\brief Parameter factor to some unit by index
 *
 * This method must be used with an advanced unit object
 */
         double getFactorToSomeUnit(std::string target,int ipar)          const {return parameter[ipar].getFactorToSomeUnit(target);}
/*!\brief Parameter factor to some unit by index, unsigned int overload
 *
 * This method must be used with an advanced unit object
 */
         double getFactorToSomeUnit(std::string target,unsigned int ipar) const {return parameter[ipar].getFactorToSomeUnit(target);}
/*!\brief Parameter factor to some unit by name
 *
 * This method must be used with an advanced unit object
 */
         double getFactorToSomeUnit(std::string target,std::string namePar)    const {return getParameter(namePar).getFactorToSomeUnit(target);}
/*!\brief Parameter change to some unit by index
 *
 * This method must be used with an advanced unit object
 */
         void changeParameterToSomeUnit(std::string target,int ipar)            {parameter[ipar].changeToSomeUnit(target);}
/*!\brief Parameter change to some unit by index, unsigned int overload.
 *
 * This method must be used with an advanced unit object.
 */
         void changeParameterToSomeUnit(std::string target,unsigned int ipar)   {parameter[ipar].changeToSomeUnit(target);}
/*!\brief Parameter change to some unit by name
 *
 * This method must be used with an advanced unit object.
 */
         void changeParameterToSomeUnit(std::string target,std::string namePar);
/*!\brief Internal unit densification. 
 *
 * All the parameters of the SimpleChem are tested against each other.
 */
         void densifyParametersUnit();
/*!\brief External unit densification. 
 *
 * All the parameters of the SimpleChem are tested against the parameters 
 * of another SimpleChem object.
 */
         void densifyParametersUnitExternal(const SimpleChem &out);
/*!\brief External unit densification.
 *
 * In this function, we densify a specific parameter with an external ParameterPhy.
 */
         bool densifyParameterToThis(int ipar, const ParameterPhy &rhs)    {return parameter[ipar].densifyUnitObject(rhs);}  
/*!\brief External unit densification.
 *
 * In this function, we densify all the parameters of
 * the SimpleChem to an external ParameterPhy.
 */
         void densifyParametersToThis(const ParameterPhy &rhs);

/*!\brief Parameters unit getter*/
         const std::vector<std::string> getParametersUnit() const;
/*!\brief Parameters unit setter*/
         void setParametersUnit(const std::vector<std::string> &units);
/*!\brief Parameters unit adder*/
         void addParametersByUnit(const std::vector<std::string> &units);

/*!\brief Parameter uncertainty pointer getter by index.*/
         CoreUnc *getUncObjectParameterPtr(int ipar)                           {return parameter[ipar].getUncertaintyObject();}
/*!\brief Parameter uncertainty pointer getter by index, unsigned int overload.*/
         CoreUnc *getUncObjectParameterPtr(unsigned int ipar)                  {return parameter[ipar].getUncertaintyObject();}
/*!\brief Parameter uncertainty pointer getter by name.*/
         CoreUnc *getUncObjectParameterPtr(std::string namePar);
/*!\brief Parameter uncertainty pointer setter by index.*/
         void setUncertaintyObjectParameter(int ipar,CoreUnc *target)          {parameter[ipar].setUncertaintyObject(target);}
/*!\brief Parameter uncertainty pointer setter by index, unsigned int overload.*/
         void setUncertaintyObjectParameter(unsigned int ipar,CoreUnc *target) {parameter[ipar].setUncertaintyObject(target);}
/*!\brief Parameter uncertainty pointer setter by name.*/
         void setUncertaintyObjectParameter(std::string namePar,CoreUnc *target);

/*!\brief Parameter uncertainty type getter by index.*/
         const std::string getParameterUncertaintyType(int ipar)               const {return parameter[ipar].typeDPar();}
/*!\brief Parameter uncertainty type getter by index, unsigned int overload.*/
         const std::string getParameterUncertaintyType(unsigned int ipar)      const {return parameter[ipar].typeDPar();}
/*!\brief Parameter uncertainty type getter by name.*/
         const std::string getParameterUncertaintyType(std::string namePar)         const {return getParameter(namePar).typeDPar();}
/*!\brief Parameter uncertainty type setter by index.*/
         void setParameterUncertaintyType(int ipar, std::string parUncType)          {parameter[ipar].setTypeDPar(parUncType);}
/*!\brief Parameter uncertainty type getter by index, unsigned int overload.*/
         void setParameterUncertaintyType(unsigned int ipar, std::string parUncType) {parameter[ipar].setTypeDPar(parUncType);}
/*!\brief Parameter uncertainty type getter by name.*/
         void setParameterUncertaintyType(std::string namePar, std::string parUncType);
/*!\brief Parameter adder with only the uncertainty type.*/
         void addParameterByUncertaintyType(std::string parUncType);

/*!\brief Parameters uncertainty type getter.*/
         const std::vector<std::string> getParametersUncertaintyType() const;
/*!\brief Parameters uncertainty type setter.*/
         void setParametersUncertaintyType(const std::vector<std::string> &parUncType);
/*!\brief Parameters uncertainty type adder.*/
         void addParametersByUncertaintyType(const std::vector<std::string> &parUncType);

/*!\brief Parameter value getter by index.*/
         double getParameterValue(int ipar, int ival)                   const   {return parameter[ipar].value(ival);}
/*!\brief Parameter value getter by index, unsigned int overload.*/
         double getParameterValue(unsigned int ipar, unsigned int ival) const   {return parameter[ipar].value(ival);}
/*!\brief Parameter value getter by name.*/
         double getParameterValue(std::string namePar, int ival)             const   {return getParameter(namePar).value(ival);}
/*!\brief Parameter value getter by name, unsigned int overload.*/
         double getParameterValue(std::string namePar, unsigned int ival)    const   {return getParameter(namePar).value(ival);}
/*!\brief Parameter value setter by index.*/
         void setParameterValue(double value,int ipar, int ival)                      {parameter[ipar].setValue(value,ival);}
/*!\brief Parameter value setter by index, unsigned int overload.*/
         void setParameterValue(double value,unsigned int ipar, unsigned int ival)    {parameter[ipar].setValue(value,ival);}
/*!\brief Parameter value setter by name.*/
         void setParameterValue(double value,std::string namePar, int ival);
/*!\brief Parameter value adder by index.*/
         void addParameterValue(double value,int ipar)                                {parameter[ipar].addValue(value);}
/*!\brief Parameter value adder by index, unsigned int overload.*/
         void addParameterValue(double value,unsigned int ipar)                       {parameter[ipar].addValue(value);}
/*!\brief Parameter value adder by name.*/
         void addParameterValue(double value,std::string namePar);

/*!\brief Parameter values getter by index.*/
         const std::vector<double> getParameterValues(int ipar)              const   {return parameter[ipar].values();}
/*!\brief Parameter values getter by index, unsigned int overload.*/
         const std::vector<double> getParameterValues(unsigned int ipar)     const   {return parameter[ipar].values();}
/*!\brief Parameter values getter by name.*/
         const std::vector<double> getParameterValues(std::string namePar)        const   {return getParameter(namePar).values();}
/*!\brief Parameter values adder by index.*/
         void addParameterValues(const std::vector<double> &vals, int ipar)          {parameter[ipar].addValues(vals);}
/*!\brief Parameter values adder by index, unsigned int overload.*/
         void addParameterValues(const std::vector<double> &vals, unsigned int ipar) {parameter[ipar].addValues(vals);}
/*!\brief Parameter values adder by name.*/
         void addParameterValues(const std::vector<double> &vals, std::string namePar);
/*\brief Parameter duplicater*/
         void duplicateParameterValue(int ival, int ntimes, int ipar)                            {parameter[ipar].duplicateValue(ival,ntimes);}
/*\brief Parameter duplicater, unsigned int overload*/
         void duplicateParameterValue(unsigned int ival, unsigned int ntimes, unsigned int ipar) {parameter[ipar].duplicateValue(ival,ntimes);}
/*\brief Parameter duplicater, by name*/
         void duplicateParameterValue(int ival, int ntimes, const std::string & namePar)              {getParameter(namePar).duplicateValue(ival,ntimes);}

/*!\brief Parameter number of values getter by index.*/
         int getParameterNValues(int ipar)                        const   {return parameter[ipar].nValues();}
/*!\brief Parameter number of values getter by index, unsigned int overload.*/
         int getParameterNValues(unsigned int ipar)               const   {return parameter[ipar].nValues();}
/*!\brief Parameter number of values getter by name.*/
         int getParameterNValues(std::string namePar)             const   {return getParameter(namePar).nValues();}
/*!\brief Parameter clearer by index.
 *
 * To clear a parameter does not erase it, it only erase
 * its content (values & name). To erase completely a parameter, use erase method.
 */
         void clearParameter(int ipar)                                          {parameter[ipar].clear();}
/*!\brief Parameter clearer by index, unsigned int overload.*/
         void clearParameter(unsigned int ipar)                                 {parameter[ipar].clear();}
/*!\brief Parameter clearer by name.*/
         void clearParameter(std::string namePar);
/*!\brief Clear all the parameters*/
         void clearParameters();
/*!\brief Parameter values clearer by index.
 *
 * This method clears only the parameter's values
 */
         void clearParameterValues(int ipar)                                    {parameter[ipar].clearValues();}
/*!\brief Parameter values clearer by index, unsigned int overload.*/
         void clearParameterValues(unsigned int ipar)                           {parameter[ipar].clearValues();}
/*!\brief Parameter values clearer by name.*/
         void clearParameterValues(std::string namePar);
/*!\brief Clear all the parameters*/
         void clearParametersValues();
/*!\brief Parameter eraser by index. Beware of mixing it with clear()*/
         void eraseParameter(int ipar)                                          {parameter.erase(parameter.begin() + ipar - 1);}
/*!\brief Parameter eraser by index, unsigned int overload.*/
         void eraseParameter(unsigned int ipar)                                 {eraseParameter((int)ipar);}
/*!\brief Parameter eraser by name.*/
         void eraseParameter(const std::string &namePar);
/*!\brief Erase all the parameters*/
         void eraseParameters()                                                 {parameter.clear();}

/*!\brief Uncertainty value of a parameter getter by index.*/
         double getParameterDValue(int ipar,int ival)                   const   {return parameter[ipar].Dvalue(ival);}
/*!\brief Uncertainty value of a parameter getter by index, unsigned int overload.*/
         double getParameterDValue(unsigned int ipar,unsigned int ival) const   {return parameter[ipar].Dvalue(ival);}
/*!\brief Uncertainty value of a parameter getter by name.*/
         double getParameterDValue(std::string namePar,int ival)             const   {return getParameter(namePar).Dvalue(ival);}
/*!\brief Uncertainty value of a parameter setter by index.*/
         void setParameterDValue(double value,int ipar,int ival)                      {parameter[ipar].setDValue(value,ival);}
/*!\brief Uncertainty value of a parameter setter by index, unsigned int overload.*/
         void setParameterDValue(double value,unsigned int ipar,unsigned int ival)    {parameter[ipar].setDValue(value,ival);}
/*!\brief Uncertainty value of a parameter setter by name.*/
         void setParameterDValue(double value,std::string namePar,int ival);
/*!\brief Uncertainty value of a parameter adder by index.*/
         void addParameterDValue(double value,int ipar)                               {parameter[ipar].addDValue(value);}
/*!\brief Uncertainty value of a parameter adder by index, unsigned int overload.*/
         void addParameterDValue(double value,unsigned int ipar)                      {parameter[ipar].addDValue(value);}
/*!\brief Uncertainty value of a parameter adder by name.*/
         void addParameterDValue(double value,std::string namePar);

/*!\brief Parameter uncertainty values getter by index.*/
         const std::vector<double> getParameterDValues(int ipar)              const   {return parameter[ipar].DValues();}
/*!\brief Parameter uncertainty values getter by index, unsigned int overload.*/
         const std::vector<double> getParameterDValues(unsigned int ipar)     const   {return parameter[ipar].DValues();}
/*!\brief Parameter uncertainty values getter by name.*/
         const std::vector<double> getParameterDValues(std::string namePar)        const   {return getParameter(namePar).DValues();}
/*!\brief Parameter uncertainty values adder by index.*/
         void addParameterDValues(const std::vector<double> &vals, int ipar)          {parameter[ipar].addDValues(vals);}
/*!\brief Parameter uncertainty values adder by index, unsigned int overload.*/
         void addParameterDValues(const std::vector<double> &vals, unsigned int ipar) {parameter[ipar].addDValues(vals);}
/*!\brief Parameter uncertainty values adder by name.*/
         void addParameterDValues(const std::vector<double> &vals, std::string namePar);
/*!\brief Parameter number of uncertainty values getter by index.*/
         int getParameterNDValues(int ipar)                        const   {return parameter[ipar].nDValues();}
/*!\brief Parameter number of uncertainty values getter by index, unsigned int overload.*/
         int getParameterNDValues(unsigned int ipar)               const   {return parameter[ipar].nDValues();}
/*!\brief Parameter number of uncertainty values getter by name.*/
         int getParameterNDValues(std::string namePar)                  const   {return getParameter(namePar).nDValues();}

/*!\brief Parameter uncertainty absolute value getter by index.*/
         double getParameterDValueAbsolute(int ipar, int ival)                  const {return parameter[ipar].getDValueAbsolute(ival);}
/*!\brief Parameter uncertainty absolute value getter by index, unsigned int overload.*/
         double getParameterDValueAbsolute(unsigned int ipar,unsigned int ival) const {return parameter[ipar].getDValueAbsolute(ival);}
/*!\brief Parameter uncertainty absolute value getter by name.*/
         double getParameterDValueAbsolute(std::string namePar, int ival)            const {return getParameter(namePar).getDValueAbsolute(ival);}
/*!\brief Parameter uncertainty absolute value getter by name, unsigned int overload.*/
         double getParameterDValueAbsolute(std::string namePar, unsigned int ival)   const {return getParameter(namePar).getDValueAbsolute(ival);}


/*\brief Parameter inserter*/
         void insertParameter(int ind, ParameterPhy &addedPar)                        {parameter.insert(parameter.begin() + ind, addedPar);}
/*\brief Parameter inserter, unsigned int overload*/
         void insertParameter(unsigned int ind, ParameterPhy &addedPar)               {parameter.insert(parameter.begin() + ind, addedPar);}

/*!\brief Parameter resizer by index*/
         void resizePar(int newSize,int ipar)                   {parameter[ipar].resizePar(newSize);}
/*!\brief Parameter resizer by index, unsigned int overload*/
         void resizePar(unsigned int newSize,unsigned int ipar) {parameter[ipar].resizePar(newSize);}
/*!\brief Parameter resizer by name*/
         void resizePar(int newSize,std::string namePar);
/*!\brief Parameter resizer by name, unsigned int overload*/
         void resizePar(unsigned int newSize,std::string namePar)    {resizePar((int)newSize,namePar);}

/*!brief Function number getter*/
         int getNFunct()                                             const {return (int)function.size();}
/*!brief Function getter*/
         PhyFunct getFunction(int i)                                 const {return function[i];}
/*!brief Function getter by name*/
         PhyFunct getFunction(const std::string &name)               const;
/*!brief Functions getter*/
         std::vector<PhyFunct> getFunctions()                        const {return function;}
/*!brief Function setter*/
         void setFunction(const PhyFunct &func,int i = 0);
/*!brief Functions setter*/
         void setFunctions(const std::vector<PhyFunct> &funcs)       {function = funcs;}
/*!brief Function adder*/
         void addFunction(const PhyFunct &func)                      {function.push_back(func);}
/*!brief Functions adder*/
         void addFunctions(const std::vector<PhyFunct> &funcs);

/*!brief HarmFunct number getter*/
         int getNHarmFunct()                                          const {return (int)harmFunct.size();}
/*!brief HarmFunct getter*/
         PhyHarmFunct getHarmFunct(int i)                             const {return harmFunct[i];}
/*!brief HarmFunct getter by name*/
         PhyHarmFunct getHarmFunct(const std::string &name)           const;
/*!brief HarmFuncts getter*/
         std::vector<PhyHarmFunct> getHarmFuncts()                    const {return harmFunct;}
/*!brief HarmFunct setter*/
         void setHarmFunct(const PhyHarmFunct &relat,int i = 0)       {harmFunct[i] = relat;}
/*!brief HarmFuncts setter*/
         void setHarmFuncts(const std::vector<PhyHarmFunct> &relats)  {harmFunct = relats;}
/*!brief HarmFunct adder*/
         void addHarmFunct(const PhyHarmFunct &rela)                  {harmFunct.push_back(rela);}
/*!brief HarmFuncts adder*/
         void addHarmFuncts(const std::vector<PhyHarmFunct> &relas);


/*\brief Parameters getter, reference version*/
         std::vector<ParameterPhy> & getParameters()                 {return parameter;}
/*\brief Parameter getter, reference version
 *
 * This is the version used for parameter calculations,
 * this version enable the parameter to be changeable.
 */
         ParameterPhy & getParameter(int ipar)                  {return parameter[ipar];}
/*\brief Parameter getter, reference version, unsigned int version*/
         ParameterPhy & getParameter(unsigned int ipar)         {return parameter[ipar];}
/*\brief Parameter getter, reference version, by name*/
         ParameterPhy & getParameter(const std::string &namePar);


/*!\brief Sends back the index of a targeted parameter/function/multifunction.
 *
 * If target is found, it sends back its index, if not
 * found, it sends back -1.
 */
         int getIndexByName(const std::string &name, const std::string &type) const;


/*!\brief showAll(std::ostream (=cout)) method.
 *
 * This method prints a full description of a SimpleChem
 * object:
 *     - Molecules involved
 *     - Description
 *     - Number of parameters
 *     - Name of parameters
 *     - Unit of parameters
 *     - Uncertainty type of parameters
 *     - Values of parameters
 *     - Uncertainty values of parameters
 */
         virtual const void showAll(std::ostream &out = std::cout)const;

/*!\brief To complete the information of a SimpleChem object with another SimpleChem object.
 *
 * This method will use the informations of the given
 * SimpleChem object to complete the stored informations.
 * If the given object has a unit object and not the
 * current one, it is copied. With a same fashion,
 * the uncertainty values and the values are compared,
 * and if needed, added.
 */
        void complete(const SimpleChem &rhs);
/*!\brief Clear all the SimpleChem object.
 *
 * It clear the molecules and description, and erase all the parameters.
 */
        void clear();

/*\brief Adress getter of this*/
        SimpleChem * getPtr()  {return this;}

/*!\brief Assignement operator.*/
        SimpleChem & operator = (const SimpleChem &rhs);
/*!\brief Alternative to SimpleChem &operator=(const SimpleChem&)*/
        void equalizeSC(const SimpleChem &rhs);               

        protected:

/*!\brief A std::vector of std::strings for the molecules concerned.*/
         std::vector<std::string> molecules;
/*!\brief A std::vector of ParameterPhy for the parameter.*/
         std::vector<ParameterPhy> parameter;
/*!\brief A std::vector of PhyFunct for coupled parameters.*/
         std::vector<PhyFunct> function;
/*!\brief A std::vector of PhyHarmFunct for related parameters.*/
         std::vector<PhyHarmFunct> harmFunct;
/*!\brief A std::string for the description.*/
         std::string description;

};
}
#endif
