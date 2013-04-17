//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
#ifndef _KINETICS_REACTION_STRUCTURE_
#define _KINETICS_REACTION_STRUCTURE_

#include "antioch/KineticsChem.hpp"

namespace Antioch{
/* \class ReactionChem
 * \brief A class to model a chemical reaction.
 *
 * To model a chemical reaction, we split the disappearing and
 * appearing matter flux:
 *   - a SimpleChem object to model the disappearing flux, the reactive part
 *      of the reaction;
 *   - a SimpleChem object per opened channel, thus a std::vector of SimpleChem
 *     objects for the appearing flux, the produced part of the reaction.
 *
 * In this fashion, a reactive channel will be found as the combination
 * of the SimpleChem object for the reactive part and the concerned SimpleChem
 * object of the produced part.
 *
 * The idea of this class is to have the information in the following manner:
 *   - The reactants and physical parameters for the global rate constant
 *      are stored in the SimpleChem object for the reactive part:
 *         - the std::vector of std::strings molecules contains the name of the reacting molecule,
 *         - the std::vector of ParameterPhy parameter contains the different needed parameters,
 *         - the std::string description contains an identifier for the kinetics model used, see
 *           the file kinetics_models.hpp for more details.
 *   - Each reactive channel is stored in the very same fashion.
 */
class ReactionChem{
  public:
/*!\brief Default constructor*/
   ReactionChem(){}
/*!\brief Copy constructor*/
   ReactionChem(const ReactionChem &rhs);
/*!\brief Default destructor*/
   ~ReactionChem(){}

/*!\brief Global rate constant setter*/
   void setGlobalRateConstant(const KineticsChem &chem);
/*!\brief Global rate constant setter for SimpleChem part*/
   void setGlobalRateConstantNonReactivePart(const SimpleChem &chem) {globalRate.setNonReactivePart(chem);}
/*!\brief Global rate constant completer for kinetics part*/
   void setGlobalRateConstantKineticsPart(const KineticsChem &chem)  {globalRate.setKineticsOnly(chem);}
/*!\brief Global rate constant getter*/
   const KineticsChem getGlobalRate()                          const {return globalRate;}
/*!\brief Global rate constant getter, SimpleChem part*/
   const SimpleChem getGlobalRateNonReactivePart()             const {return *(static_cast<const SimpleChem*>(&globalRate));}
/*!\brief Global rate constant pointer getter*/
   KineticsChem * getPtrGlobalRate()                                 {return &globalRate;}

/*!\brief Reactive channel setter*/
   void setReactiveChannel(const SimpleChem &chem, int ichan)          {reactiveChannel[ichan].equalizeSC(chem);}
/*!\brief Reactive channel setter, unsigned int overload*/
   void setReactiveChannel(const SimpleChem &chem, unsigned int ichan) {setReactiveChannel(chem,(int)ichan);}
/*!\brief Reactive channel adder*/
   void addReactiveChannel(const SimpleChem &chem)                     {reactiveChannel.push_back(chem);}
/*!\brief Empty reactive channel adder*/
   void addEmptyReactiveChannel()                                      {reactiveChannel.push_back(SimpleChem());}
/*!\brief Reactive channel getter*/
   SimpleChem getReactiveChannel(int ichan = 0)                  const {return reactiveChannel[ichan];}
/*!\brief Reactive channel getter, unsigned int overload*/
   SimpleChem getReactiveChannel(unsigned int ichan)             const {return getReactiveChannel((int)ichan);}
/*!\brief Reactive channel pointer getter*/
   SimpleChem * getPtrReactiveChannel(int ichan = 0)                   {return reactiveChannel[ichan].getPtr();}
/*!\brief Reactive channel pointer getter, unsigned int overload*/
   SimpleChem * getPtrReactiveChannel(unsigned int ichan)              {return reactiveChannel[ichan].getPtr();}

/*!\brief Reactive part setter*/
   void setAllReactiveChannels(const std::vector<SimpleChem> &chem);
/*!\brief Reactive part getter*/
   std::vector<SimpleChem> getAllReactiveChannels()             const {return reactiveChannel;}
/*!\brief Test if the given set of molecules is a reactive channel*/
   bool GotTheseProducts(const std::vector<std::string> &mols)  const;
/*!\brief Index of reactive channel defined by the set of molecules, returns -1 if not found*/
   int WhereTheseProducts(const std::vector<std::string> &mols)  const;

/*!\brief Set the ReactionChem object*/
   void setReaction(const ReactionChem &reac);


/*!\brief Number of branching ratios getter*/
   int getNChannels()                                          const {return reactiveChannel.size();}

/*!\brief Total number of parameters getters
 *
 * This method returns the number of parameters of the global
 * rate constant plus the number of parameters of each
 * branching ratio.
 */
   int nPars()              const;

////////////////////////////////////////////////
////Kinetics rate constant part
///////////////////////////////////////////////

/*!\brief Total number of kinetics parameters getters
 *
 * Only the kinetics part is concerned here.
 *
 */
   int nKinPars()                                                   const {return globalRate.nKinPars();}

/*!\brief Kinetics reaction setter*/
    void setGlobalRateConstantKinetics(Process * kin)                    {globalRate.setKinetics(kin);}
/*!\brief Chemical process of rate constant getter.*/
   const std::string getGlobalRateConstantKineticsProcess()         const {return globalRate.getKineticsProcess();}
/*!\brief Chemical process of rate constant setter.*/
   void setGlobalRateConstantKineticsProcess(const std::string &kinProc)  {globalRate.setKineticsProcess(kinProc);}
/*!\brief Kinetics model of rate constant getter.*/
   const std::string getGlobalRateConstantKineticsModel()           const {return globalRate.getKineticsModel();}
/*!\brief Kinetics model of rate constant setter.*/
   void setGlobalRateConstantKineticsModel(const std::string &kinModel)   {globalRate.setKineticsModel(kinModel);}

/*!\brief Set the values of temperature*/
   void setKineticsTemperature(ParameterPhy * temp)          {globalRate.setKineticsTemperature(temp);}
/*!\brief Set the values of concentration*/
   void setKineticsConcentration(ParameterPhy * conc)        {globalRate.setKineticsConcentration(conc);}
/*!\brief Rate constant calculations, all is internally defined*/
   std::vector<ParameterPhy> allRateConstant()         const {return globalRate.allRateConstant();}
/*!\brief Rate constant calculations, all is internally defined*/
   ParameterPhy getRateConstant(int nk = 0)            const {return globalRate.rateConstant(nk);}
/*!\brief Rate constant calculations, temp is given*/
   double getRateConstantT(double t, int i = 0)        const {return globalRate.getRateConstantT(t,i);}
/*!\brief Rate constant calculations, for falloff
 *
 * Here the concentration is given, the temperature
 * must be defined. The concentration will be considered
 * in default unit: mol/cm3
 */
    ParameterPhy getRateConstantM(double m, int i = 0)            {return globalRate.rateConstantM(m, i);}
/*!\brief Rate constant calculations, for falloff
 *
 * Here all is given, the temperature and
 * the concentration will be considered
 * in default unit, respectively in K and mol/cm3.
 */
    double getRateConstantTM(double t, double m, int i = 0) const {return globalRate.getRateConstantTM(t, m, i);}
/*!\brief Rate constant calculations, for falloff
 *
 * Here the temperature is given, the concentration
 * must be defined. Useful for k(M) calculation at
 * fixed temperature (in K, always).
 */
    ParameterPhy rateConstantT(double t)                          {return globalRate.rateConstantT(t);}
/*!\brief Sets the kinetics part
 *
 * The parameters are stored within the kinetics object,
 * sorted and then the kinetics is initialized.
 */
    void setKinetics(const std::string &kinProc, const std::string &kinModel, const std::vector<ParameterPhy> &pars);

    void showKinetics(std::ostream &out = std::cout)  const {globalRate.showKinetics(out);}

/*!\brief Gives a statistical report over the value of a partial rate constant
 *\todo redo this method, its outdated
 *
 * This method returns in a std::string the statistical summary of
 * a channel. It is given in the form of a tabular:
 *
 *  parameter name | parameter unit | 5% quantile | 95% quantile | 50% quantile | mean | standard deviation
 * 
 * If the parameter if the preexponential factor, then its value is multiplied
 * by the value of the considered branching ratio.
 */
   const std::string report(int ichan);

///////////////////////////////////////
///////////////////////////////////////

/*!\brief Gives the statistical report over all the reactive channels.*/
   const std::string reportAll();
/*!\brief Scan all the parameters to densify to the given ParameterPhy.*/
   void densifyParametersToThis(const ParameterPhy &par);
/*!\brief Scan all the parameters to densify them between each other, within SimpleChem object.*/
   void densifyAllUnitsIntra();
/*!\brief Scan all the parameters to densify them between each other, wherever they are.*/
   void densifyAllUnitsExtra();
/*!\brief Scan all the parameters to densify them to the given ReactionChem.*/
   void densifyAllUnitsExternal(ReactionChem &out);

// Global rate constant
/*!\brief showAll() of the global rate constant.*/
   void GlobalRateConstantShowAll(std::ostream &out = std::cout) const;
/*!\brief clear() of the global rate constant.*/
   void GlobalRateConstantClear()                               {globalRate.clear();}
/*!\brief Number of reactants getter.*/
   int NReactants()                                       const {return globalRate.getNMolecules();}
/*!\brief Name of reactants getter.*/
   const std::string Reactant(int ireac)                  const {return globalRate.getMolecule(ireac);}
/*!\brief Name of reactants getter, unsigned int overload.*/
   const std::string Reactant(unsigned int ireac)         const {return globalRate.getMolecule(ireac);}
/*!\brief Reactant adder.*/
   void addReactant(std::string reac)                           {globalRate.addMolecule(reac);}
/*!\brief Reactant setter.*/
   void setReactant(std::string reac,int ireac)                 {globalRate.setMolecule(reac,ireac);}
/*!\brief Reactant setter, unsigned int overload.*/
   void setReactant(std::string reac,unsigned int ireac)        {globalRate.setMolecule(reac,ireac);}

/*!\brief Reactants getter.*/
   const std::vector<std::string> Reactants()             const {return globalRate.getMolecules();}
/*!\brief Reactants setter.*/
   void setReactants(const std::vector<std::string> &reacs)     {globalRate.setMolecules(reacs);}
/*!\brief Reactants tester.*/
   bool areSameReactants(const std::vector<std::string> &mols) const {return globalRate.areSameMolecules(mols);}
/*!\brief Reactants tester, returns the index of the reactant if found, -1 if not found.*/
   int isInReactants(const std::string &mol) const {return globalRate.IsInMolecules(mol);}

/*!\brief Number of parameter of global rate constant getter.*/
   int GlobalRateConstantNParameters()                    const {return globalRate.getNPar();}

/*!\brief Kinetics model getter for a reactive channel.*/
   const std::string GlobalRateConstantDescription()                       const {return globalRate.getDescription();}
// parameters of rate constant
/*!\brief Add an empty parameter to the global rate constant.*/
   void addGlobalRateConstantEmptyParameter()                                    {globalRate.addEmptyParameter();}
/*!\brief Index finder of parameter.*/
   int GlobalRateConstantGetIndexByName(const std::string &name, const std::string &type) const {return globalRate.getIndexByName(name,type);}
/*!\brief Parameter adder to the global rate constant.*/
   void addGlobalRateConstantParameter(const ParameterPhy &par)                  {globalRate.addParameter(par);}
/*!\brief Parameter getter.*/
   const ParameterPhy getGlobalRateConstantParameter(int ipar)             const {return globalRate.getParameter(ipar);}
/*!\brief Parameter getter, unsigned int overload.*/
   const ParameterPhy getGlobalRateConstantParameter(unsigned int ipar)    const {return globalRate.getParameter(ipar);}
/*!\brief Parameter getter by name.*/
   const ParameterPhy getGlobalRateConstantParameter(std::string namePar)  const {return globalRate.getParameter(namePar);}
/*!\brief Name getter.*/
   const std::string GlobalRateConstantParName(int i)           const {return globalRate.getParameterName(i);}
/*!\brief Name getter, unsigned int overload.*/
   const std::string GlobalRateConstantParName(unsigned int i)  const {return globalRate.getParameterName(i);}
/*!\brief Names getter.*/
   const std::vector<std::string> GlobalRateConstantParsName()       const {return globalRate.getParametersName();}
/*!\brief Name setter.*/
   void setGlobalRateConstantParName(std::string name,int i)          {globalRate.setParameterName(name,i);}
/*!\brief Name setter, unsigned int overload.*/
   void setGlobalRateConstantParName(std::string name,unsigned int i) {globalRate.setParameterName(name,i);}

/*!\brief Parameters getter, reference version for ParameterPhy calculation*/
   std::vector<ParameterPhy> &getGlobalRateConstantParameters()       {return globalRate.getParameters();}
/*!\brief Parameter getter, reference version for ParameterPhy calculation*/
   ParameterPhy &getGlobalRateConstantParameter(int ipar)             {return globalRate.getParameter(ipar);}
/*!\brief Parameter getter, reference version for ParameterPhy calculation, unsigned int overload*/
   ParameterPhy &getGlobalRateConstantParameter(unsigned int ipar)    {return globalRate.getParameter(ipar);}
/*!\brief Parameter getter, reference version for ParameterPhy calculation, by name*/
   ParameterPhy &getGlobalRateConstantParameter(const std::string &namePar) {return globalRate.getParameter(namePar);}

// units
/*!\brief Unit object pointer of global rate constant getter.*/
   Units *getGlobalRateConstantUnitObjectParameterPtr(int ipar)                       {return globalRate.getUnitObjectParameterPtr(ipar);}
/*!\brief Unit object pointer of global rate constant getter, unsigned int overload.*/
   Units *getGlobalRateConstantUnitObjectParameterPtr(unsigned int ipar)              {return globalRate.getUnitObjectParameterPtr(ipar);}
/*!\brief Unit object pointer of global rate constant getter by name.*/
   Units *getGlobalRateConstantUnitObjectParameterPtr(std::string namePar)            {return globalRate.getUnitObjectParameterPtr(namePar);}
/*!\brief Unit object pointer of global rate constant setter.*/
   void setGlobalRateConstantUnitObjectParameter(int ipar, Units *target)             {globalRate.setUnitObjectParameter(ipar,target);}
/*!\brief Unit object pointer of global rate constant setter, unsigned int overload.*/
   void setGlobalRateConstantUnitObjectParameter(unsigned int ipar,Units *target)     {globalRate.setUnitObjectParameter(ipar,target);}
/*!\brief Unit object pointer of global rate constant setter by name.*/
   void setGlobalRateConstantUnitObjectParameter(std::string namePar,Units *target)   {globalRate.setUnitObjectParameter(namePar,target);}

/*!\brief Global rate constant densifier to a ParameterPhy.*/
   bool densifyGlobalRateConstantParameterToThis(int ipar,const ParameterPhy &par) {return globalRate.densifyParameterToThis(ipar,par);}

/*!\brief Unit of parameter getter for global rate constant.*/
   const std::string getGlobalRateConstantParUnit(int i)              const {return globalRate.getParameterUnit(i);}
/*!\brief Unit of parameter getter for global rate constant, unsigned int overload.*/
   const std::string getGlobalRateConstantParUnit(unsigned int i)     const {return globalRate.getParameterUnit(i);}
/*!\brief Unit of parameter getter by name for global rate constant.*/
   const std::string getGlobalRateConstantParUnit(std::string namePar)     const {return globalRate.getParameterUnit(namePar);}
/*!\brief Units of parameters getter for global rate constant.*/
   const std::vector<std::string> getGlobalRateConstantParsUnit()          const {return globalRate.getParametersUnit();}
/*!\brief SI factor getter for global rate constant.*/
   double getGlobalRateConstantParSIFactor(int i)          const {return globalRate.getFactorToSI(i);}
/*!\brief SI factor getter for global rate constant, unsigned int overload.*/
   double getGlobalRateConstantParSIFactor(unsigned int i) const {return globalRate.getFactorToSI(i);}
/*!\brief SI factor getter by name for global rate constant.*/
   double getGlobalRateConstantParSIFactor(std::string namePar) const {return globalRate.getFactorToSI(namePar);}
/*!\brief Factor to any unit getter for global rate constant.*/
   double getGlobalRateConstantParFactorToSomeUnit(std::string target,int ipar)          const {return globalRate.getFactorToSomeUnit(target,ipar);}
/*!\brief Factor to any unit getter for global rate constant, unsigned int overload.*/
   double getGlobalRateConstantParFactorToSomeUnit(std::string target,unsigned int ipar) const {return globalRate.getFactorToSomeUnit(target,ipar);}
/*!\brief Factor to any unit getter by name for global rate constant.*/
   double getGlobalRateConstantParFactorToSomeUnit(std::string target,std::string namePar)    const {return globalRate.getFactorToSomeUnit(target,namePar);}
/*!\brief Unit changer for global rate constant.*/
   void   changeGlobalRateConstantParameterToSomeUnit(std::string target,int ipar)            {globalRate.changeParameterToSomeUnit(target,ipar);}
/*!\brief Unit changer for global rate constant, unsigned int overload.*/
   void   changeGlobalRateConstantParameterToSomeUnit(std::string target,unsigned int ipar)   {globalRate.changeParameterToSomeUnit(target,ipar);}
/*!\brief Unit changer by name for global rate constant.*/
   void   changeGlobalRateConstantParameterToSomeUnit(std::string target,std::string namePar)      {globalRate.changeParameterToSomeUnit(target,namePar);}

/*!\brief Unit setter for global rate constant.*/
   void setGlobalRateConstantParUnit(std::string TargetUnit,int i)          {globalRate.setParameterUnit(i,TargetUnit);}
/*!\brief Unit setter for global rate constant, unsigned int overload.*/
   void setGlobalRateConstantParUnit(std::string TargetUnit,unsigned int i) {globalRate.setParameterUnit(i,TargetUnit);}
/*!\brief Unit setter by name for global rate constant.*/
   void setGlobalRateConstantParUnit(std::string TargetUnit,std::string namePar) {globalRate.setParameterUnit(namePar,TargetUnit);}

/*!\brief Uncertainty type getter for global rate constant.*/
   const std::string GlobalRateConstantParUncertaintyType(int ipar)          const {return globalRate.getParameterUncertaintyType(ipar);}
/*!\brief Uncertainty type getter for global rate constant, unsigned int overload.*/
   const std::string GlobalRateConstantParUncertaintyType(unsigned int ipar) const {return globalRate.getParameterUncertaintyType(ipar);}
/*!\brief Uncertainty type getter by name for global rate constant.*/
   const std::string GlobalRateConstantParUncertaintyType(std::string namePar)    const {return globalRate.getParameterUncertaintyType(namePar);}
/*!\brief Uncertainty type getter of all parameters for global rate constant.*/
   const std::vector<std::string> GlobalRateConstantParsUncertaintyType()         const {return globalRate.getParametersUncertaintyType();}
/*!\brief Uncertainty type setter for global rate constant.*/
   void setGlobalRateConstantParUncertaintyType(std::string uncType, int i)           {globalRate.setParameterUncertaintyType(i,uncType);}
/*!\brief Uncertainty type setter for global rate constant, unsigned int overload.*/
   void setGlobalRateConstantParUncertaintyType(std::string uncType, unsigned int i)  {globalRate.setParameterUncertaintyType(i,uncType);}
/*!\brief Uncertainty type setter by name for global rate constant.*/
   void setGlobalRateConstantParUncertaintyType(std::string uncType, std::string namePar)  {globalRate.setParameterUncertaintyType(namePar,uncType);}

/*!\brief Parameter value getter for global rate constant.*/
   double GlobalRateConstantParValue(int ipar, int ival)                     const {return globalRate.getParameterValue(ipar,ival);}
/*!\brief Parameter value getter by name for global rate constant.*/
   double GlobalRateConstantParValue(std::string namePar, int ival)               const {return globalRate.getParameterValue(namePar,ival);}
/*!\brief Parameter value getter by name for global rate constant, unsigned int overload.*/
   double GlobalRateConstantParValue(std::string namePar, unsigned int ival)      const {return globalRate.getParameterValue(namePar,ival);}
/*!\brief Parameter value setter for global rate constant.*/
   void setGlobalRateConstantParValue(double val, int ipar, int ival)                    {globalRate.setParameterValue(val,ipar,ival);}
/*!\brief Parameter value setter for global rate constant, unsigned int overload.*/
   void setGlobalRateConstantParValue(double val, unsigned int ipar, unsigned int ival)  {globalRate.setParameterValue(val,ipar,ival);}
/*!\brief Parameter value setter by name for global rate constant.*/
   void setGlobalRateConstantParValue(double val, std::string namePar, int ival)              {globalRate.setParameterValue(val,namePar,ival);}
/*!\brief Parameter value adder for global rate constant.*/
   void addGlobalRateConstantParValue(double val, int ipar)                              {globalRate.addParameterValue(val,ipar);}
/*!\brief Parameter value adder for global rate constant, unsigned int overload.*/
   void addGlobalRateConstantParValue(double val, unsigned int ipar)                     {globalRate.addParameterValue(val,ipar);}
/*!\brief Parameter value adder by name for global rate constant.*/
   void addGlobalRateConstantParValue(double val, std::string namePar)                        {globalRate.addParameterValue(val,namePar);}

/*!\brief Values getter for global rate constant.*/
   const std::vector<double> GlobalRateConstantParValues(int ipar)                   const {return globalRate.getParameterValues(ipar);}
/*!\brief Values getter for global rate constant, unsigned int overload.*/
   const std::vector<double> GlobalRateConstantParValues(unsigned int ipar)          const {return globalRate.getParameterValues(ipar);}
/*!\brief Values getter by name for global rate constant.*/
   const std::vector<double> GlobalRateConstantParValues(std::string namePar)             const {return globalRate.getParameterValues(namePar);}
/*!\brief Values adder for global rate constant.*/
   void addGlobalRateConstantParValues(const std::vector<double> &pars, int ipar)          {globalRate.addParameterValues(pars,ipar);}
/*!\brief Values adder for global rate constant, unsigned int overload.*/
   void addGlobalRateConstantParValues(const std::vector<double> &pars, unsigned int ipar) {globalRate.addParameterValues(pars,ipar);}
/*!\brief Values adder by name for global rate constant.*/
   void addGlobalRateConstantParValues(const std::vector<double> &pars, std::string namePar)    {globalRate.addParameterValues(pars,namePar);}
/*!\brief Values duplicater for global rate constant.*/
   void duplicateGlobalRateConstantParValue(int ival, int ntimes, int ipar)                            {globalRate.duplicateParameterValue(ival,ntimes,ipar);}
/*!\brief Values duplicater for global rate constant., unsigned int overload*/
   void duplicateGlobalRateConstantParValue(unsigned int ival, unsigned int ntimes, unsigned int ipar) {globalRate.duplicateParameterValue(ival,ntimes,ipar);}
/*!\brief Values duplicater by name for global rate constant.*/
   void duplicateGlobalRateConstantParValue(int ival, int ntimes, std::string namePar)                      {globalRate.duplicateParameterValue(ival,ntimes,namePar);}


/*!\brief Number of values getter for global rate constant.*/
   int GlobalRateConstantParNValues(int ipar = 0)                               const {return globalRate.getParameterNValues(ipar);}
/*!\brief Number of values getter for global rate constant, unsigned int overload.*/
   int GlobalRateConstantParNValues(unsigned int ipar)                          const {return globalRate.getParameterNValues(ipar);}
/*!\brief Number of values getter by name for global rate constant.*/
   int GlobalRateConstantParNValues(std::string namePar)                             const {return globalRate.getParameterNValues(namePar);}
/*!\brief Parameter clearer for global rate constant.*/
   void GlobalRateConstantClearParameter(int ipar)                                    {globalRate.clearParameter(ipar);}
/*!\brief Parameter clearer for global rate constant, unsigned int overload.*/
   void GlobalRateConstantClearParameter(unsigned int ipar)                           {globalRate.clearParameter(ipar);}
/*!\brief Parameter clearer by name for global rate constant.*/
   void GlobalRateConstantClearParameter(std::string namePar)                              {globalRate.clearParameter(namePar);}
/*!\brief All parameters clearer for global rate constant.*/
   void GlobalRateConstantClearParameters()                                           {globalRate.clearParameters();}
/*!\brief Parameter values clearer for global rate constant.*/
   void GlobalRateConstantClearParameterValues(int ipar)                              {globalRate.clearParameterValues(ipar);}
/*!\brief Parameter values clearer for global rate constant, unsigned int overload.*/
   void GlobalRateConstantClearParameterValues(unsigned int ipar)                     {globalRate.clearParameterValues(ipar);}
/*!\brief Parameter clearer by name for global rate constant.*/
   void GlobalRateConstantClearParameterValues(std::string namePar)                        {globalRate.clearParameterValues(namePar);}
/*!\brief All parameters clearer for global rate constant.*/
   void GlobalRateConstantClearParametersValues()                                     {globalRate.clearParametersValues();}

/*!\brief Uncertainty value getter for global rate constant.*/
   double GlobalRateConstantParDValue(int ipar, int ival)                    const {return globalRate.getParameterDValue(ipar,ival);}
/*!\brief Uncertainty value getter for global rate constant, unsigned int overload.*/
   double GlobalRateConstantParDValue(unsigned int ipar, unsigned int ival)  const {return globalRate.getParameterDValue(ipar,ival);}
/*!\brief Uncertainty value getter by name for global rate constant.*/
   double GlobalRateConstantParDValue(std::string namePar, int ival)              const {return globalRate.getParameterDValue(namePar,ival);}
/*!\brief Uncertainty value setter for global rate constant.*/
   void setGlobalRateConstantParDValue(double val, int ipar, int ival)                   {globalRate.setParameterDValue(val,ipar,ival);}
/*!\brief Uncertainty value setter for global rate constant, unsigned overload.*/
   void setGlobalRateConstantParDValue(double val, unsigned int ipar, unsigned int ival) {globalRate.setParameterDValue(val,ipar,ival);}
/*!\brief Uncertainty value setter by name for global rate constant.*/
   void setGlobalRateConstantParDValue(double val, std::string namePar, int ival)             {globalRate.setParameterDValue(val,namePar,ival);}
/*!\brief Uncertainty value adder for global rate constant.*/
   void addGlobalRateConstantParDValue(double val, int ipar)                             {globalRate.addParameterDValue(val,ipar);}
/*!\brief Uncertainty value adder for global rate constant, unsigned int overload.*/
   void addGlobalRateConstantParDValue(double val, unsigned int ipar)                    {globalRate.addParameterDValue(val,ipar);}
/*!\brief Uncertainty value adder by name for global rate constant.*/
   void addGlobalRateConstantParDValue(double val, std::string namePar)                       {globalRate.addParameterDValue(val,namePar);}

/*!\brief Uncertainty values getter for global rate constant.*/
   const std::vector<double> GlobalRateConstantParDValues(int ipar)                   const {return globalRate.getParameterDValues(ipar);}
/*!\brief Uncertainty values getter for global rate constant, unsigned int overload.*/
   const std::vector<double> GlobalRateConstantParDValues(unsigned int ipar)          const {return globalRate.getParameterDValues(ipar);}
/*!\brief Uncertainty values getter by name for global rate constant.*/
   const std::vector<double> GlobalRateConstantParDValues(std::string namePar)             const {return globalRate.getParameterDValues(namePar);}
/*!\brief Uncertainty values adder for global rate constant.*/
   void addGlobalRateConstantParDValues(const std::vector<double> &pars, int ipar)          {globalRate.addParameterDValues(pars,ipar);}
/*!\brief Uncertainty values adder for global rate constant, unsigned int overload.*/
   void addGlobalRateConstantParDValues(const std::vector<double> &pars, unsigned int ipar) {globalRate.addParameterDValues(pars,ipar);}
/*!\brief Uncertainty values adder by name for global rate constant.*/
   void addGlobalRateConstantParDValues(const std::vector<double> &pars, std::string namePar)    {globalRate.addParameterDValues(pars,namePar);}
/*!\brief Number of uncertainty values getter for global rate constant.*/
   int GlobalRateConstantParNDValues(int ipar)                                   const {return globalRate.getParameterNDValues(ipar);}
/*!\brief Number of uncertainty values getter for global rate constant, unsigned int overload.*/
   int GlobalRateConstantParNDValues(unsigned int ipar)                          const {return globalRate.getParameterNDValues(ipar);}
/*!\brief Number of uncertainty values getter by name for global rate constant.*/
   int GlobalRateConstantParNDValues(std::string namePar)                             const {return globalRate.getParameterNDValues(namePar);}

/*!\brief Absolute value of uncertainty getter for global rate constant.*/
   double GlobalRateConstantParDValueAbsolute(int ipar, int ival)                    const {return globalRate.getParameterDValueAbsolute(ipar,ival);}
/*!\brief Absolute value of uncertainty getter for global rate constant, unsigned int overload.*/
   double GlobalRateConstantParDValueAbsolute(unsigned int ipar, unsigned int ival)  const {return globalRate.getParameterDValueAbsolute(ipar,ival);}
/*!\brief Absolute value of uncertainty getter by name for global rate constant.*/
   double GlobalRateConstantParDValueAbsolute(std::string namePar, int ival)              const {return globalRate.getParameterDValueAbsolute(namePar,ival);}
/*!\brief Absolute value of uncertainty getter by name for global rate constant, unsigned int overload.*/
   double GlobalRateConstantParDValueAbsolute(std::string namePar, unsigned int ival)     const {return globalRate.getParameterDValueAbsolute(namePar,ival);}

/*!\brief Parameter resizer for global rate constant.*/
   void GlobalRateConstantResize(int newSize,int ipar)                                     {globalRate.resizePar(newSize,ipar);}
/*!\brief Parameter resizer for global rate constant, unsigned int overload.*/
   void GlobalRateConstantResize(unsigned int newSize,unsigned int ipar)                   {globalRate.resizePar(newSize,ipar);}
/*!\brief Parameter resizer by name for global rate constant.*/
   void GlobalRateConstantResize(unsigned int newSize,std::string namePar)                      {globalRate.resizePar(newSize,namePar);}

/*!\brief Parameter inserter.*/
   void GlobalRateConstantInsertParameter(int ind, ParameterPhy &newPar)                   {globalRate.insertParameter(ind,newPar);}
/*!\brief Parameter inserter, unsigned int overload.*/
   void GlobalRateConstantInsertParameter(unsigned int ind, ParameterPhy &newPar)          {globalRate.insertParameter(ind,newPar);}

// Reactive channels
/*!\brief Reactive channel showAll().*/
   void ReactiveChannelShowAll(int ichan = 0, std::ostream &out = std::cout)      const {out << "Reactive channel " << ichan;reactiveChannel[ichan].showAll(out);}
/*!\brief Reactive channel showAll(), unsigned int overload.*/
   void ReactiveChannelShowAll(unsigned int ichan, std::ostream &out = std::cout) const {out << "Reactive channel " << ichan;reactiveChannel[ichan].showAll(out);}
/*!\brief Clear of the reactive channel.
 *
 * Clearing the reactive channel does not erase it. It
 * makes it empty.
 */
   void ReactiveChannelClear(int ichan = 0)              {reactiveChannel[ichan].clear();}
/*!\brief Clear of the reactive channel, unsigned int overload.*/
   void ReactiveChannelClear(unsigned int ichan)         {reactiveChannel[ichan].clear();}
/*!\brief Clear all reactive channels.*/
   void ReactiveChannelsClearAll();
/*!\brief Erase the reactive channel ichan.*/
   void ReactiveChannelErase(int ichan = 0)              {reactiveChannel.erase(reactiveChannel.begin() + ichan - 1);}
/*!\brief Erase reactive channel ichan, unsigned int overload.*/
   void ReactiveChannelErase(unsigned int ichan)         {ReactiveChannelErase((int)ichan);}
/*!\brief Erase all reactive channels.*/
   void ReactiveChannelsEraseAll()                       {reactiveChannel.clear();}

/*!\brief Number of products getter.*/
   int NProducts(int ichan = 0)                                const {return reactiveChannel[ichan].getNMolecules();}
/*!\brief Number of products getter, unsigned int overload.*/
   int NProducts(unsigned int ichan)                           const {return reactiveChannel[ichan].getNMolecules();}
/*!\brief Product name getter.*/
   const std::string Product(int iprod,int ichan = 0)               const {return reactiveChannel[ichan].getMolecule(iprod);}
/*!\brief Product name getter, unsigned int overload.*/
   const std::string Product(unsigned int iprod,unsigned int ichan) const {return reactiveChannel[ichan].getMolecule(iprod);}
/*!\brief Product adder.*/
   void addProduct(std::string prod,int ichan = 0)                        {reactiveChannel[ichan].addMolecule(prod);}
/*!\brief Product adder, unsigned int overload.*/
   void addProduct(std::string prod,unsigned int ichan)                   {reactiveChannel[ichan].addMolecule(prod);}
/*!\brief Product setter.*/
   void setProduct(std::string prod,int iprod,int ichan = 0)              {reactiveChannel[ichan].setMolecule(prod,iprod);}
/*!\brief Product setter, unsigned int overload.*/
   void setProduct(std::string prod,unsigned int iprod,unsigned int ichan){reactiveChannel[ichan].setMolecule(prod,iprod);}

/*!\brief Products getter.*/
   const std::vector<std::string> Products(int ichan = 0)                         const {return reactiveChannel[ichan].getMolecules();}
/*!\brief Products getter, unsigned int overload.*/
   const std::vector<std::string> Products(unsigned int ichan)                    const {return reactiveChannel[ichan].getMolecules();}
/*!\brief Products setter.*/
   void setProducts(const std::vector<std::string> &prods, int ichan = 0)               {reactiveChannel[ichan].setMolecules(prods);}
/*!\brief Products setter, unsigned int overload.*/
   void setProducts(const std::vector<std::string> &prods, unsigned int ichan)          {reactiveChannel[ichan].setMolecules(prods);}
/*!\brief Products tester.*/
   bool areSameProducts(const std::vector<std::string> &mols, int ichan = 0)      const {return reactiveChannel[ichan].areSameMolecules(mols);}
/*!\brief Products tester, unsigned int overload.*/
   bool areSameProducts(const std::vector<std::string> &mols, unsigned int ichan) const {return reactiveChannel[ichan].areSameMolecules(mols);}
/*!\brief Products tester, returns the index of the product of found, -1 if not.*/
   int isInProducts(const std::string &mol, int ichan = 0)                   const {return reactiveChannel[ichan].IsInMolecules(mol);}
/*!\brief Products tester, returns the index of the product of found, -1 if not, unsigned int overload.*/
   int isInProducts(const std::string &mol, unsigned int ichan)                   const {return reactiveChannel[ichan].IsInMolecules(mol);}

/*!\brief Kinetics model getter for a reactive channel.*/
   const std::string ReactiveChannelKineticsType(int ichan = 0)                           const {return reactiveChannel[ichan].getDescription();}
/*!\brief Kinetics model getter for a reactive channel, unsigned int overload.*/
   const std::string ReactiveChannelKineticsType(unsigned int ichan)                      const {return reactiveChannel[ichan].getDescription();}
/*!\brief Kinetics model setter for a reactive channel.*/
   void setReactiveChannelKineticsType(std::string kinModel, int ichan = 0)                     {reactiveChannel[ichan].setDescription(kinModel);}
/*!\brief Kinetics model setter for a reactive channel, unsigned int overload.*/
   void setReactiveChannelKineticsType(std::string kinModel, unsigned int ichan)                {reactiveChannel[ichan].setDescription(kinModel);}

/*!\brief Branching ratios getter*/
   const std::vector<ParameterPhy> ReactiveChannelBranchingRatios();
/*!\brief Branching ratio getter*/
   ParameterPhy ReactiveChannelBranchingRatio(int ichan = 0);
/*!\brief Branching ratio ParameterPhy  for a reactive channel, reference overload for ParameterPhy calculations.*/
//   ParameterPhy &ReactiveChannelBranchingRatio(int ichan = 0);
/*!\brief Value of branching ratio getter for a reactive channel.*/
   double ReactiveChannelBranchingRatioValue(int ival, int ichan = 0);
/*!\brief Value of branching ratio getter for a reactive channel, unsigned int overload.*/
   double ReactiveChannelBranchingRatioValue(unsigned int ival, unsigned  int ichan)       {return ReactiveChannelBranchingRatioValue((int)ival,(int)ichan);};
/*!\brief Number of parameters getter for a reactive channel.*/
   int ReactiveChannelNParameters(int ichan = 0)                                     const {return reactiveChannel[ichan].getNPar();}
/*!\brief Number of parameters getter for a reactive channel, unsigned int overload.*/
   int ReactiveChannelNParameters(unsigned int ichan)                                const {return reactiveChannel[ichan].getNPar();}

/*!\brief Empty parameter adder for a reactive channel.*/
   void addReactiveChannelEmptyParameter(int ichan = 0)                                {reactiveChannel[ichan].addEmptyParameter();}
/*!\brief Empty parameter adder for a reactive channel, unsigned int overload.*/
   void addReactiveChannelEmptyParameter(unsigned int ichan)                           {reactiveChannel[ichan].addEmptyParameter();}
/*!\brief Parameter index getter for a reactive channel.*/
   int ReactiveChannelIndexByName(std::string name,const std::string &type, int ichan = 0) const {return reactiveChannel[ichan].getIndexByName(name,type);}
/*!\brief Parameter index getter for a reactive channel, unsigned int overload.*/
   int ReactiveChannelIndexParameterByName(std::string name,const std::string &type, unsigned int ichan) const {return reactiveChannel[ichan].getIndexByName(name,type);}
/*!\brief Parameter adder for a reactive channel.*/
   void ReactiveChannelAddParameter(const ParameterPhy &par,int ichan = 0)                    {reactiveChannel[ichan].addParameter(par);}
/*!\brief Parameter adder for a reactive channel, unsigned int overload.*/
   void ReactiveChannelAddParameter(const ParameterPhy &par,unsigned int ichan)               {reactiveChannel[ichan].addParameter(par);}
/*!\brief Parameter getter for a reactive channel.*/
   const ParameterPhy getReactiveChannelParameter(int ipar,int ichan = 0)               const {return reactiveChannel[ichan].getParameter(ipar);}
/*!\brief Parameter getter for a reactive channel, unsigned int overload.*/
   const ParameterPhy getReactiveChannelParameter(unsigned int ipar,unsigned int ichan) const {return reactiveChannel[ichan].getParameter(ipar);}
/*!\brief Parameter getter by name for a reactive channel.*/
   const ParameterPhy getReactiveChannelParameter(std::string namePar, int ichan = 0)        const {return reactiveChannel[ichan].getParameter(namePar);}
/*!\brief Parameter getter by name for a reactive channel, unsigned int overload.*/
   const ParameterPhy getReactiveChannelParameter(std::string namePar, unsigned int ichan)   const {return reactiveChannel[ichan].getParameter(namePar);}
/*!\brief Parameter getter for a reactive channel, reference for ParameterPhy calculations * .*/
   ParameterPhy &getReactiveChannelParameter(int ipar,int ichan = 0)                    {return reactiveChannel[ichan].getParameter(ipar);}
/*!\brief Parameter getter for a reactive channel, unsigned int overload, reference for ParameterPhy calculations.*/
   ParameterPhy &getReactiveChannelParameter(unsigned int ipar,unsigned int ichan)      {return reactiveChannel[ichan].getParameter(ipar);}
/*!\brief Parameter getter by name for a reactive channel, reference for ParameterPhy calculations.*/
   ParameterPhy &getReactiveChannelParameter(std::string namePar, int ichan = 0)        {return reactiveChannel[ichan].getParameter(namePar);}
/*!\brief Parameter getter by name for a reactive channel, unsigned int overload, reference for ParameterPhy calculations.*/
   ParameterPhy &getReactiveChannelParameter(std::string namePar, unsigned int ichan)   {return reactiveChannel[ichan].getParameter(namePar);}
/*!\brief Parameter name getter for a reactive channel.*/
   const std::string ReactiveChannelParName(int i, int ichan = 0)                            const {return reactiveChannel[ichan].getParameterName(i);}
/*!\brief Parameter name getter for a reactive channel, unsigned int overload.*/
   const std::string ReactiveChannelParName(unsigned int i, unsigned int ichan)              const {return reactiveChannel[ichan].getParameterName(i);}
/*!\brief Parameters names getter for a reactive channel.*/
   const std::vector<std::string> ReactiveChannelParsName(int ichan = 0)                          const {return reactiveChannel[ichan].getParametersName();}
/*!\brief Parameters names getter for a reactive channel, unsigned int overload.*/
   const std::vector<std::string> ReactiveChannelParsName(unsigned int ichan)                     const {return reactiveChannel[ichan].getParametersName();}
/*!\brief Parameter name setter for a reactive channel.*/
   void setReactiveChannelParName(std::string name,int i, int ichan = 0)                {reactiveChannel[ichan].setParameterName(name,i);}
/*!\brief Parameter name setter for a reactive channel, unsigned int overload.*/
   void setReactiveChannelParName(std::string name,unsigned int i, unsigned int ichan)  {reactiveChannel[ichan].setParameterName(name,i);}

// units
/*!\brief Unit pointer getter for a reactive channel.*/
   Units *getReactiveChannelUnitObjectParameterPtr(int ipar,int ichan = 0)                        {return reactiveChannel[ichan].getUnitObjectParameterPtr(ipar);}
/*!\brief Unit pointer getter for a reactive channel, unsigned int overload.*/
   Units *getReactiveChannelUnitObjectParameterPtr(unsigned int ipar,unsigned int ichan)          {return reactiveChannel[ichan].getUnitObjectParameterPtr(ipar);}
/*!\brief Unit pointer getter by name for a reactive channel.*/
   Units *getReactiveChannelUnitObjectParameterPtr(std::string namePar,int ichan = 0)                  {return reactiveChannel[ichan].getUnitObjectParameterPtr(namePar);}
/*!\brief Unit pointer getter by name for a reactive channel, unsigned int overload.*/
   Units *getReactiveChannelUnitObjectParameterPtr(std::string namePar,unsigned int ichan)             {return reactiveChannel[ichan].getUnitObjectParameterPtr(namePar);}
/*!\brief Unit pointer setter for a reactive channel.*/
   void setReactiveChannelUnitObjectParameter(int ipar,Units *target,int ichan = 0)               {reactiveChannel[ichan].setUnitObjectParameter(ipar,target);}
/*!\brief Unit pointer setter for a reactive channel, unsigned int overload.*/
   void setReactiveChannelUnitObjectParameter(unsigned int ipar,Units *target,unsigned int ichan) {reactiveChannel[ichan].setUnitObjectParameter(ipar,target);}
/*!\brief Unit pointer setter by name for a reactive channel.*/
   void setReactiveChannelUnitObjectParameter(std::string namePar,Units *target,int ichan = 0)         {reactiveChannel[ichan].setUnitObjectParameter(namePar,target);}
/*!\brief Unit pointer setter by name for a reactive channel, unsigned int overload.*/
   void setReactiveChannelUnitObjectParameter(std::string namePar,Units *target,unsigned int ichan)    {reactiveChannel[ichan].setUnitObjectParameter(namePar,target);}

/*!\brief Parameter densifier to a given ParameterPhy for a reactive channel.*/
   bool densifyReactiveChannelParameterToThis(int ipar,const ParameterPhy &par, int ichan = 0)       {return reactiveChannel[ichan].densifyParameterToThis(ipar,par);}
/*!\brief Parameter densifier to a given ParameterPhy for a reactive channel, unsigned int overload.*/
   bool densifyReactiveChannelParameterToThis(int ipar,const ParameterPhy &par, unsigned int ichan)  {return reactiveChannel[ichan].densifyParameterToThis(ipar,par);}

/*!\brief Parameter unit getter for a reactive channel.*/
   const std::string ReactiveChannelParUnit(int ipar, int ichan = 0)               const {return reactiveChannel[ichan].getParameterUnit(ipar);}
/*!\brief Parameter unit getter for a reactive channel, unsigned int overload.*/
   const std::string ReactiveChannelParUnit(unsigned int ipar, unsigned int ichan) const {return reactiveChannel[ichan].getParameterUnit(ipar);}
/*!\brief Parameter unit getter by name for a reactive channel.*/
   const std::string ReactiveChannelParUnit(std::string namePar, int ichan = 0)         const {return reactiveChannel[ichan].getParameterUnit(namePar);}
/*!\brief Parameter unit getter by name for a reactive channel, unsigned int overload.*/
   const std::string ReactiveChannelParUnit(std::string namePar, unsigned int ichan)    const {return reactiveChannel[ichan].getParameterUnit(namePar);}
/*!\brief Parameters units getter for a reactive channel.*/
   const std::vector<std::string> ReactiveChannelParsUnit(int ichan = 0)                const {return reactiveChannel[ichan].getParametersUnit();}
/*!\brief Parameters units getter for a reactive channel, unsigned int overload.*/
   const std::vector<std::string> ReactiveChannelParsUnit(unsigned int ichan)           const {return reactiveChannel[ichan].getParametersUnit();}
/*!\brief Factor to SI unit getter for a reactive channel.*/
   double ReactiveChannelParSIFactor(int ipar,int ichan = 0)                  const {return reactiveChannel[ichan].getFactorToSI(ipar);}
/*!\brief Factor to SI unit getter for a reactive channel, unsigned int overload.*/
   double ReactiveChannelParSIFactor(unsigned int ipar,unsigned int ichan)    const {return reactiveChannel[ichan].getFactorToSI(ipar);}
/*!\brief Factor to SI unit getter by name for a reactive channel.*/
   double ReactiveChannelParSIFactor(std::string namePar, int ichan = 0)           const {return reactiveChannel[ichan].getFactorToSI(namePar);}
/*!\brief Factor to SI unit getter by name for a reactive channel, unsigned int overload.*/
   double ReactiveChannelParSIFactor(std::string namePar, unsigned int ichan)      const {return reactiveChannel[ichan].getFactorToSI(namePar);}

/*!\brief Factor to any unit getter for a reactive channel.*/
   double ReactiveChannelParFactorToSomeUnit(std::string target,int ipar,int ichan = 0)               const {return reactiveChannel[ichan].getFactorToSomeUnit(target,ipar);}
/*!\brief Factor to any unit getter for a reactive channel, unsigned int overload.*/
   double ReactiveChannelParFactorToSomeUnit(std::string target,unsigned int ipar, unsigned int ichan)const {return reactiveChannel[ichan].getFactorToSomeUnit(target,ipar);}
/*!\brief Factor to any unit getter by name for a reactive channel.*/
   double ReactiveChannelParFactorToSomeUnit(std::string target,std::string namePar, int ichan = 0)     const {return reactiveChannel[ichan].getFactorToSomeUnit(target,namePar);}
/*!\brief Factor to any unit getter by name for a reactive channel, unsigned int overload.*/
   double ReactiveChannelParFactorToSomeUnit(std::string target,std::string namePar, unsigned int ichan)const {return reactiveChannel[ichan].getFactorToSomeUnit(target,namePar);}
/*!\brief Factor to any unit getter for a reactive channel.*/
   void ReactiveChannelChangeParameterToSomeUnit(std::string target,int ipar,int ichan = 0)              {reactiveChannel[ichan].changeParameterToSomeUnit(target,ipar);}
/*!\brief Factor to any unit getter for a reactive channel, unsigned int overload.*/
   void ReactiveChannelChangeParameterToSomeUnit(std::string target,unsigned int ipar,unsigned int ichan){reactiveChannel[ichan].changeParameterToSomeUnit(target,ipar);}
/*!\brief Factor to any unit getter by name for a reactive channel.*/
   void ReactiveChannelChangeParameterToSomeUnit(std::string target,std::string namePar,int ichan = 0)        {reactiveChannel[ichan].changeParameterToSomeUnit(target,namePar);}
/*!\brief Factor to any unit getter by name for a reactive channel, unsigned int overload.*/
   void ReactiveChannelChangeParameterToSomeUnit(std::string target,std::string namePar,unsigned int ichan)   {reactiveChannel[ichan].changeParameterToSomeUnit(target,namePar);}

/*!\brief Unit setter for a reactive channel.*/
   void setReactiveChannelParUnit(int ipar, std::string TargetUnit, int ichan = 0)               {reactiveChannel[ichan].setParameterUnit(ipar,TargetUnit);}
/*!\brief Unit setter for a reactive channel, unsigned int overload.*/
   void setReactiveChannelParUnit(unsigned int ipar, std::string TargetUnit, unsigned int ichan) {reactiveChannel[ichan].setParameterUnit(ipar,TargetUnit);}
/*!\brief Unit setter by name for a reactive channel.*/
   void setReactiveChannelParUnit(std::string namePar, std::string TargetUnit, int ichan = 0)         {reactiveChannel[ichan].setParameterUnit(namePar,TargetUnit);}
/*!\brief Unit setter by name for a reactive channel, unsigned int overload.*/
   void setReactiveChannelParUnit(std::string namePar, std::string TargetUnit, unsigned int ichan)    {reactiveChannel[ichan].setParameterUnit(namePar,TargetUnit);}

/*!\brief Parameter uncertainty type getter for a reactive channel.*/
   const std::string ReactiveChannelParUncertaintyType(int ipar, int ichan = 0)               const {return reactiveChannel[ichan].getParameterUncertaintyType(ipar);}
/*!\brief Parameter uncertainty type getter for a reactive channel, unsigned int overload.*/
   const std::string ReactiveChannelParUncertaintyType(unsigned int ipar, unsigned int ichan) const {return reactiveChannel[ichan].getParameterUncertaintyType(ipar);}
/*!\brief Parameter uncertainty type getter by name for a reactive channel.*/
   const std::string ReactiveChannelParUncertaintyType(std::string namePar, int ichan = 0)         const {return reactiveChannel[ichan].getParameterUncertaintyType(namePar);}
/*!\brief Parameter uncertainty type getter by name for a reactive channel, unsigned int overload.*/
   const std::string ReactiveChannelParUncertaintyType(std::string namePar, unsigned int ichan)    const {return reactiveChannel[ichan].getParameterUncertaintyType(namePar);}
/*!\brief Parameters uncertainty types getter for a reactive channel.*/
   const std::vector<std::string> ReactiveChannelParsUncertaintyType(int ichan = 0)                const {return reactiveChannel[ichan].getParametersUncertaintyType();}
/*!\brief Parameters uncertainty types getter for a reactive channel, unsigned int overload.*/
   const std::vector<std::string> ReactiveChannelParsUncertaintyType(unsigned int ichan)           const {return reactiveChannel[ichan].getParametersUncertaintyType();}
/*!\brief Parameter uncertainty type setter for a reactive channel.*/
   void setReactiveChannelParUncertaintyType(std::string uncType,int ipar, int ichan = 0)
                                                                {reactiveChannel[ichan].setParameterUncertaintyType(ipar,uncType);}
/*!\brief Parameter uncertainty type setter for a reactive channel, unsigned int overload.*/
   void setReactiveChannelParUncertaintyType(std::string uncType,unsigned int ipar, unsigned int ichan)
                                                                {reactiveChannel[ichan].setParameterUncertaintyType(ipar,uncType);}
/*!\brief Parameter uncertainty type setter by name for a reactive channel.*/
   void setReactiveChannelParUncertaintyType(std::string uncType,std::string namePar, int ichan = 0)
                                                                {reactiveChannel[ichan].setParameterUncertaintyType(namePar,uncType);}
/*!\brief Parameter uncertainty type setter by name for a reactive channel, unsigned int overload.*/
   void setReactiveChannelParUncertaintyType(std::string uncType,std::string namePar, unsigned int ichan)
                                                                {reactiveChannel[ichan].setParameterUncertaintyType(namePar,uncType);}

/*!\brief Parameter value getter for a reactive channel.*/
   double ReactiveChannelParValue(int ipar, int ival, int ichan = 0)                         const {return reactiveChannel[ichan].getParameterValue(ipar,ival);}
/*!\brief Parameter value getter for a reactive channel, unsigned int overload.*/
   double ReactiveChannelParValue(unsigned int ipar, unsigned int ival, unsigned int ichan)  const {return reactiveChannel[ichan].getParameterValue(ipar,ival);}
/*!\brief Parameter value getter by name for a reactive channel.*/
   double ReactiveChannelParValue(std::string namePar, int ival, int ichan = 0)                   const {return reactiveChannel[ichan].getParameterValue(namePar,ival);}
/*!\brief Parameter value getter by name for a reactive channel, unsigned int overload.*/
   double ReactiveChannelParValue(std::string namePar, int ival, unsigned int ichan)              const {return reactiveChannel[ichan].getParameterValue(namePar,ival);}

/*!\brief Parameter value setter for a reactive channel.*/
   void setReactiveChannelParValue(double val, int ipar, int ival, int ichan = 0)                        {reactiveChannel[ichan].setParameterValue(val,ipar,ival);}
/*!\brief Parameter value setter for a reactive channel, unsigned int overload.*/
   void setReactiveChannelParValue(double val, unsigned int ipar, unsigned int ival, unsigned int ichan) {reactiveChannel[ichan].setParameterValue(val,ipar,ival);}
/*!\brief Parameter value setter by name for a reactive channel.*/
   void setReactiveChannelParValue(double val, std::string namePar, int ival, int ichan = 0)                  {reactiveChannel[ichan].setParameterValue(val,namePar,ival);}
/*!\brief Parameter value setter by name for a reactive channel, unsigned int overload.*/
   void setReactiveChannelParValue(double val, std::string namePar, int ival, unsigned int ichan)             {reactiveChannel[ichan].setParameterValue(val,namePar,ival);}
/*!\brief Parameter value setter by name for a reactive channel, unsigned int overload.*/
   void setReactiveChannelParValue(double val, std::string namePar, unsigned int ival, unsigned int ichan)    {reactiveChannel[ichan].setParameterValue(val,namePar,ival);}
/*!\brief Parameter value adder for a reactive channel.*/
   void addReactiveChannelParValue(double val, int ipar, unsigned int ichan)                             {reactiveChannel[ichan].addParameterValue(val,ipar);}
/*!\brief Parameter value adder for a reactive channel, unsigned int overload.*/
   void addReactiveChannelParValue(double val, unsigned int ipar, unsigned int ichan)                    {reactiveChannel[ichan].addParameterValue(val,ipar);}
/*!\brief Parameter value adder by name for a reactive channel.*/
   void addReactiveChannelParValue(double val, std::string namePar, int ichan)                                {reactiveChannel[ichan].addParameterValue(val,namePar);}
/*!\brief Parameter value adder by name for a reactive channel, unsigned int overload.*/
   void addReactiveChannelParValue(double val, std::string namePar, unsigned int ichan)                       {reactiveChannel[ichan].addParameterValue(val,namePar);}

/*!\brief Parameter value duplicater for a reactive channel.*/
   void duplicateReactiveChannelParValue(int ival, int ntimes, int ipar, unsigned int ichan)                             {reactiveChannel[ichan].duplicateParameterValue(ival,ntimes,ipar);}
/*!\brief Parameter value duplicater for a reactive channel, unsigned int overload.*/
   void duplicateReactiveChannelParValue(unsigned int ival, unsigned int ntimes, unsigned int ipar, unsigned int ichan)  {reactiveChannel[ichan].duplicateParameterValue(ival,ntimes,ipar);}
/*!\brief Parameter value duplicater by name for a reactive channel.*/
   void duplicateReactiveChannelParValue(int ival, int ntimes, std::string namePar, int ichan)                                {reactiveChannel[ichan].duplicateParameterValue(ival,ntimes,namePar);}
/*!\brief Parameter value duplicater by name for a reactive channel, unsigned int overload.*/
   void duplicateReactiveChannelParValue(int ival, int ntimes, std::string namePar, unsigned int ichan)                       {reactiveChannel[ichan].duplicateParameterValue(ival,ntimes,namePar);}

/*!\brief Parameter values getter for a reactive channel.*/
   const std::vector<double> ReactiveChannelParValues(int ipar, int ichan = 0)                  const {return reactiveChannel[ichan].getParameterValues(ipar);}
/*!\brief Parameter values getter for a reactive channel, unsigned int overload.*/
   const std::vector<double> ReactiveChannelParValues(unsigned int ipar, unsigned int ichan)    const {return reactiveChannel[ichan].getParameterValues(ipar);}
/*!\brief Parameter values getter by name for a reactive channel.*/
   const std::vector<double> ReactiveChannelParValues(std::string namePar, int ichan = 0)            const {return reactiveChannel[ichan].getParameterValues(namePar);}
/*!\brief Parameter values getter by name for a reactive channel, unsigned int overload.*/
   const std::vector<double> ReactiveChannelParValues(std::string namePar, unsigned int ichan)       const {return reactiveChannel[ichan].getParameterValues(namePar);}

/*!\brief Parameter values adder for a reactive channel.*/
   void addReactiveChannelParValues(const std::vector<double> &pars, int ipar, int ichan = 0)               {reactiveChannel[ichan].addParameterValues(pars,ipar);}
/*!\brief Parameter values adder for a reactive channel, unsigned int overload.*/
   void addReactiveChannelParValues(const std::vector<double> &pars, unsigned int ipar, unsigned int ichan) {reactiveChannel[ichan].addParameterValues(pars,ipar);}
/*!\brief Parameter values adder by name for a reactive channel.*/
   void addReactiveChannelParValues(const std::vector<double> &pars, std::string namePar, int ichan = 0)         {reactiveChannel[ichan].addParameterValues(pars,namePar);}
/*!\brief Parameter values adder by name for a reactive channel, unsigned int overload.*/
   void addReactiveChannelParValues(const std::vector<double> &pars, std::string namePar, unsigned int ichan)    {reactiveChannel[ichan].addParameterValues(pars,namePar);}

/*!\brief Number of values getter for a reactive channel.*/
   int ReactiveChannelParNValues(int ipar, int ichan = 0)                            const {return reactiveChannel[ichan].getParameterNValues(ipar);}
/*!\brief Number of values getter for a reactive channel, unsigned int overload.*/
   int ReactiveChannelParNValues(unsigned int ipar, unsigned int ichan)              const {return reactiveChannel[ichan].getParameterNValues(ipar);}
/*!\brief Number of values getter by name for a reactive channel.*/
   int ReactiveChannelParNValues(std::string namePar, int ichan = 0)                      const {return reactiveChannel[ichan].getParameterNValues(namePar);}
/*!\brief Number of values getter by name for a reactive channel, unsigned int overload.*/
   int ReactiveChannelParNValues(std::string namePar, unsigned int ichan)                 const {return reactiveChannel[ichan].getParameterNValues(namePar);}
/*!\brief Parameter clearer for a reactive channel.*/
   void ReactiveChannelClearParameter(int ipar,int ichan = 0)                              {reactiveChannel[ichan].clearParameter(ipar);}
/*!\brief Parameter clearer for a reactive channel, unsigned int overload.*/
   void ReactiveChannelClearParameter(unsigned int ipar,unsigned int ichan)                {reactiveChannel[ichan].clearParameter(ipar);}
/*!\brief Parameter clearer by name for a reactive channel.*/
   void ReactiveChannelClearParameter(std::string namePar,int ichan = 0)                        {reactiveChannel[ichan].clearParameter(namePar);}
/*!\brief Parameter clearer by name for a reactive channel, unsigned int overload.*/
   void ReactiveChannelClearParameter(std::string namePar,unsigned int ichan)                   {reactiveChannel[ichan].clearParameter(namePar);}
/*!\brief All parameters clearer for a reactive channel.*/
   void ReactiveChannelClearParameters(int ichan = 0)                                      {reactiveChannel[ichan].clearParameters();}
/*!\brief All parameters clearer for a reactive channel, unsigned int overload.*/
   void ReactiveChannelClearParameters(unsigned int ichan)                                 {reactiveChannel[ichan].clearParameters();}
/*!\brief All parameters clearer for all reactive channels.*/
   void ReactiveChannelsClearParameters();
/*!\brief Parameter values clearer for a reactive channel.*/
   void ReactiveChannelClearParameterValues(int ipar,int ichan = 0)                              {reactiveChannel[ichan].clearParameterValues(ipar);}
/*!\brief Parameter values clearer for a reactive channel, unsigned int overload.*/
   void ReactiveChannelClearParameterValues(unsigned int ipar,unsigned int ichan)                {reactiveChannel[ichan].clearParameterValues(ipar);}
/*!\brief Parameter values clearer by name for a reactive channel.*/
   void ReactiveChannelClearParameterValues(std::string namePar,int ichan = 0)                        {reactiveChannel[ichan].clearParameterValues(namePar);}
/*!\brief Parameter values clearer by name for a reactive channel, unsigned int overload.*/
   void ReactiveChannelClearParameterValues(std::string namePar,unsigned int ichan)                   {reactiveChannel[ichan].clearParameterValues(namePar);}
/*!\brief All parameters values clearer for a reactive channel.*/
   void ReactiveChannelClearParametersValues(int ichan = 0)                                      {reactiveChannel[ichan].clearParametersValues();}
/*!\brief All parameters values clearer for a reactive channel, unsigned int overload.*/
   void ReactiveChannelClearParametersValues(unsigned int ichan)                                 {reactiveChannel[ichan].clearParametersValues();}
/*!\brief All parameters values clearer for all reactive channels.*/
   void ReactiveChannelsClearParametersValues();

/*!\brief Parameter uncertainty value getter for a reactive channel.*/
   double ReactiveChannelParDValue(int ipar, int ival, int ichan = 0) const 
                                                        {return reactiveChannel[ichan].getParameterDValue(ipar,ival);}
/*!\brief Parameter uncertainty value getter for a reactive channel, unsigned int overload.*/
   double ReactiveChannelParDValue(unsigned int ipar, unsigned int ival, unsigned int ichan) const 
                                                        {return reactiveChannel[ichan].getParameterDValue(ipar,ival);}
/*!\brief Parameter uncertainty value getter by name for a reactive channel.*/
   double ReactiveChannelParDValue(std::string namePar, int ival, int ichan = 0) const 
                                                        {return reactiveChannel[ichan].getParameterDValue(namePar,ival);}
/*!\brief Parameter uncertainty value getter by name for a reactive channel, unsigned int overload.*/
   double ReactiveChannelParDValue(std::string namePar, int ival, unsigned int ichan) const 
                                                        {return reactiveChannel[ichan].getParameterDValue(namePar,ival);}

/*!\brief Parameter uncertainty value setter for a reactive channel.*/
   void setReactiveChannelParDValue(double val, int ipar, int ival, int ichan = 0)                        {reactiveChannel[ichan].setParameterDValue(val,ipar,ival);}
/*!\brief Parameter uncertainty value setter for a reactive channel, unsigned int overload.*/
   void setReactiveChannelParDValue(double val, unsigned int ipar, unsigned int ival, unsigned int ichan) {reactiveChannel[ichan].setParameterDValue(val,ipar,ival);}
/*!\brief Parameter uncertainty value setter by name for a reactive channel.*/
   void setReactiveChannelParDValue(double val, std::string namePar, int ival, int ichan = 0)                  {reactiveChannel[ichan].setParameterDValue(val,namePar,ival);}
/*!\brief Parameter uncertainty value setter by name for a reactive channel, unsigned int overload.*/
   void setReactiveChannelParDValue(double val, std::string namePar, int ival, unsigned int ichan)             {reactiveChannel[ichan].setParameterDValue(val,namePar,ival);}
/*!\brief Parameter uncertainty value adder for a reactive channel.*/
   void addReactiveChannelParDValue(double val, int ipar, int ichan = 0)                                  {reactiveChannel[ichan].addParameterDValue(val,ipar);}
/*!\brief Parameter uncertainty value adder for a reactive channel, unsigned int overload.*/
   void addReactiveChannelParDValue(double val, unsigned int ipar, unsigned int ichan)                    {reactiveChannel[ichan].addParameterDValue(val,ipar);}
/*!\brief Parameter uncertainty value adder by name for a reactive channel.*/
   void addReactiveChannelParDValue(double val, std::string namePar, int ichan = 0)                            {reactiveChannel[ichan].addParameterDValue(val,namePar);}
/*!\brief Parameter uncertainty value adder by name for a reactive channel, unsigned int overload.*/
   void addReactiveChannelParDValue(double val, std::string namePar, unsigned int ichan)                       {reactiveChannel[ichan].addParameterDValue(val,namePar);}

/*!\brief Parameter uncertainty values getter for a reactive channel.*/
   const std::vector<double> ReactiveChannelParDValues(int ipar, int ichan = 0)                  const {return reactiveChannel[ichan].getParameterDValues(ipar);}
/*!\brief Parameter uncertainty values getter for a reactive channel, unsigned int overload.*/
   const std::vector<double> ReactiveChannelParDValues(unsigned int ipar, unsigned int ichan)    const {return reactiveChannel[ichan].getParameterDValues(ipar);}
/*!\brief Parameter uncertainty values getter by name for a reactive channel.*/
   const std::vector<double> ReactiveChannelParDValues(std::string namePar, int ichan = 0)            const {return reactiveChannel[ichan].getParameterDValues(namePar);}
/*!\brief Parameter uncertainty values getter by name for a reactive channel, unsigned int overload.*/
   const std::vector<double> ReactiveChannelParDValues(std::string namePar, unsigned int ichan)       const {return reactiveChannel[ichan].getParameterDValues(namePar);}

/*!\brief Parameter uncertainty values adder for a reactive channel.*/
   void addReactiveChannelParDValues(const std::vector<double> &pars, int ipar, int ichan = 0)               {reactiveChannel[ichan].addParameterDValues(pars,ipar);}
/*!\brief Parameter uncertainty values adder for a reactive channel, unsigned int overload.*/
   void addReactiveChannelParDValues(const std::vector<double> &pars, unsigned int ipar, unsigned int ichan) {reactiveChannel[ichan].addParameterDValues(pars,ipar);}
/*!\brief Parameter uncertainty values adder by name for a reactive channel.*/
   void addReactiveChannelParDValues(const std::vector<double> &pars, std::string namePar, int ichan = 0)         {reactiveChannel[ichan].addParameterDValues(pars,namePar);}
/*!\brief Parameter uncertainty values adder by name for a reactive channel, unsigned int overload.*/
   void addReactiveChannelParDValues(const std::vector<double> &pars, std::string namePar, unsigned int ichan)    {reactiveChannel[ichan].addParameterDValues(pars,namePar);}

/*!\brief Number of uncertainty values getter for a reactive channel.*/
   int ReactiveChannelParNDValues(int ipar, int ichan = 0)                            const {return reactiveChannel[ichan].getParameterNDValues(ipar);}
/*!\brief Number of uncertainty values getter for a reactive channel, unsigned int overload.*/
   int ReactiveChannelParNDValues(unsigned int ipar, unsigned int ichan)              const {return reactiveChannel[ichan].getParameterNDValues(ipar);}
/*!\brief Number of uncertainty values getter by name for a reactive channel.*/
   int ReactiveChannelParNDValues(std::string namePar, int ichan = 0)                      const {return reactiveChannel[ichan].getParameterNDValues(namePar);}
/*!\brief Number of uncertainty values getter by name for a reactive channel, unsigned int overload.*/
   int ReactiveChannelParNDValues(std::string namePar, unsigned int ichan)                 const {return reactiveChannel[ichan].getParameterNDValues(namePar);}

/*!\brief Absolute uncertainty values getter for a reactive channel.*/
   double ReactiveChannelParDValueAbsolute(int ipar, int ival,int ichan = 0)                const {return reactiveChannel[ichan].getParameterDValueAbsolute(ipar,ival);}
/*!\brief Absolute uncertainty values getter for a reactive channel, unsigned int overload.*/
   double ReactiveChannelParDValueAbsolute(unsigned int ipar, unsigned int ival, unsigned int ichan)
                                                                                             const {return reactiveChannel[ichan].getParameterDValueAbsolute(ipar,ival);}
/*!\brief Absolute uncertainty values getter by name for a reactive channel.*/
   double ReactiveChannelParDValueAbsolute(std::string namePar, int ival,int ichan = 0)      const {return reactiveChannel[ichan].getParameterDValueAbsolute(namePar,ival);}
/*!\brief Absolute uncertainty values getter by name for a reactive channel, unsigned int overload.*/
   double ReactiveChannelParDValueAbsolute(std::string namePar, int ival,unsigned int ichan) const {return reactiveChannel[ichan].getParameterDValueAbsolute(namePar,ival);}

/*!\brief Adds a reactive channel with a branching ratio of one with a none uncertainty.*/
  void addAbsoluteOneBranchingRatio(const std::vector<std::string> &mols);

/*!\brief Parameter resizer for a branching ratio.*/
  void ReactiveChannelResize(int newSize,int ipar,int ichan = 0)                          {reactiveChannel[ichan].resizePar(newSize,ipar);}
/*!\brief Parameter resizer for a branching ratio, unsigned int overload.*/
  void ReactiveChannelResize(unsigned int newSize,unsigned int ipar,unsigned int ichan)   {reactiveChannel[ichan].resizePar(newSize,ipar);}
/*!\brief Parameter resizer by name for a branching ratio.*/
  void ReactiveChannelResize(unsigned int newSize,std::string namePar,int ichan = 0)           {reactiveChannel[ichan].resizePar(newSize,namePar);}
/*!\brief Parameter resizer by name for a branching ratio, unsigned int overload.*/
  void ReactiveChannelResize(unsigned int newSize,std::string namePar,unsigned int ichan)      {reactiveChannel[ichan].resizePar(newSize,namePar);}

/*!\brief Parameter inserter.*/
   void ReactiveChannelInsertParameter(int ind, ParameterPhy &newPar, int ichan = 0)                    {reactiveChannel[ichan].insertParameter(ind,newPar);}
/*!\brief Parameter inserter, unsigned int overload.*/
   void ReactiveChannelInsertParameter(unsigned int ind, ParameterPhy &newPar, unsigned int ichan)      {reactiveChannel[ichan].insertParameter(ind,newPar);}

/*!\brief Returns the equation of the reaction.*/
   const std::string ReactionEquation() const;
/*!\brief Returns the equation of one channel of the reaction.*/
   const std::string ReactionChannelEquation(int ichan) const;
/*!\brief Prints the equation of the reaction.*/
   const void ShowReactionResume(std::ostream &out = std::cout) const {out << ReactionEquation();}
/*!\brief showAll() for the whole ReactionChem object.*/
   virtual const void showAll(std::ostream &out = std::cout) const;

/*!\brief Clears all the parameters of the ReactionChem object.
 *
 * The global rate constant is emptied, all the reactive channels
 * are erased.
 */
   void clearAll();


/*!\brief Sets all the parameter of the ReactionChem object to the same size.*/
  void resizeAll(int newSize);
/*!\brief Full resizer overload*/
  void resizeAll(unsigned int newSize)     {resizeAll((int)newSize);}

/*!\brief Id getter*/
  const std::string getID()          const {return id;}
/*!\brief Id setter*/
   void setID(const std::string &nameReac) {id = nameReac;}

  private:

/*!\brief Branching ratio calculator
 *
 * As of now (10 oct 2012) the only calculations pattern is a
 * multiplication of the ParameterPhy values. This multiplication
 * is to code for the different level of Dirichlet imbrication
 * possible.
 */
  void makeBranchingRatio(int ichan);
/*!\brief Branching ratio calculator for all reactive channels*/
  void makeBranchingRatios();

  protected:
/*!\brief A SimpleChem object for the global rate constant.*/
   KineticsChem globalRate;
/*!\brief A std::vector of SimpleChem object for the reactive channels.*/
   std::vector<SimpleChem> reactiveChannel;
/*!\brief A std::vector of ParameterPhy object for the processed branching ratios.*/
   std::vector<ParameterPhy> branchingRatios;
/*!\brief A std::string to name the reaction*/
   std::string id;
};
}
#endif
