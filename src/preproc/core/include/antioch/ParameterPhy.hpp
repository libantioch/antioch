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
#ifndef _PHYSICAL_PARAMETER_
#define _PHYSICAL_PARAMETER_

//Antioch
#include "antioch/Error.hpp"
#include "antioch/CoreUnc.hpp"
#include "antioch/Units.hpp"
#include "antioch/unit_defs.hpp"
#include "antioch/PDF.hpp"
#include "antioch/AtomicParameter.hpp"

//C++

namespace Antioch{

/*!\file PhysicalParameter.hpp
 * \brief ParameterPhy class
 * \author \me
 *
 * Contains the ParameterPhy class. 
 *
 * \class ParameterPhy
 * \todo Go through the constructors again, purge a little,
 * some are useless, and certainly adapt to user-friendly
 * assumptions:
 *   - no unit provided is ""
 *   - no uncertainty provided is CORE_UNCERTAINTY_TYPE_NONE
 *   - this implies that we should never have a NULL pointer for
 *     unit or uncertainty
 *   - What about pdf pointer?
 *
 * \brief This class groups an Unit class, an
 * Uncertainty class, a Pdf class and an AtomicParameter class
 *
 * A physical parameter is a parameter associated
 * with an uncertainty and a unit.
 * Separating the three makes it more portable, and easier
 * to deal with ParameterPhy arithmetics, as we deal of
 * each of those instances separatly.
 *
 * This class allows to combine parameters with a full
 * uncertainty and unit treatment, which depend on the
 * objects used. Theses objects are pointers to be able
 * to use any object that derives from the CoreUnit and
 * the CoreUnc objects.
 *
 * The arithmetics rely on the definition of the operators,
 * for instance, using the core objects for the unit and
 * the uncertainty:
 * \f[
 *      \varphi_1 + \varphi_2 \Leftrightarrow
 *                      \left\{\begin{array}{l}
 *                              unit_1 \mathrm{\ is\ egal\ } unit_1\\
 *                              par_1 + par_2\\
 *                              unc_1 + unc_2 \Leftrightarrow \sqrt{unc_1^2 + unc_2^2}
 *                             \end{array}\right.
 * \f]
 * with \f$\varphi\f$ being the ParameterPhy, \f$unit\f$ the unit object,
 * \f$unc\f$ the uncertainty object and \f$par\f$ the AtomicParameter
 * entity. See the class CoreUnc for the definition of the 
 * CoreUnc operator+(const CoreUnc&) const (this is a combined uncertainty,
 * not a simple sum).
 *
 * The arithmetics is described by:
 * \f[\renewcommand{\arraystretch}{1.5}
 *    \begin{array}{r@{\Leftrightarrow}l}
 *      \varphi_1 + \varphi_2 &
 *                      \left\{\begin{array}{ll}
 *                              unit_1 \mathrm{\ is\ egal\ } unit_1\\
 *                              par_1 + par_2\\
 *                              f_{add}^{uncertainty}(unc_1,unc_2)
 *                             \end{array}\right.  \tabularnewline[25pt]
 *      \varphi_1 - \varphi_2 &
 *                      \left\{\begin{array}{l}
 *                              unit_1 \mathrm{\ is\ egal\ } unit_1\\
 *                              par_1 - par_2\\
 *                              f_{sub}^{uncertainty}(unc_1,unc_2)
 *                             \end{array}\right.  \tabularnewline[25pt]
 *      \varphi_1 \times \varphi_2 &
 *                      \left\{\begin{array}{l}
 *                              unit_1 + unit_1\\
 *                              par_1 \times par_2\\
 *                              f_{mult}^{uncertainty}(unc_1,unc_2)
 *                             \end{array}\right.  \tabularnewline[25pt]
 *      \frac{\varphi_1}{\varphi_2} &
 *                      \left\{\begin{array}{l}
 *                              unit_1 - unit_1\\
 *                              \frac{par_1}{par_2}\\
 *                              f_{div}^{uncertainty}(unc_1,unc_2)
 *                             \end{array}\right.  \tabularnewline[25pt]
 *      f(\varphi_1,\dots,\varphi_n) & 
 *                      \left\{\begin{array}{l}
 *                             f_{c}^{unit}(unit_1,\dots,unit_2)\\
 *                             f(par_1,\dots,par_n)\\
 *                             f_{c}^{uncertainty}(unc_1,\dots,unc_n)
 *                             \end{array}\right.
 *    \end{array}
 *                                      
 * \f]
 *
 * The combined uncertainty calculations depends on the chosen uncertainty object.
 * The file Phymath.hpp contains the declaration of the implemented functions.
 * The functions \f$f_{add}^{uncertainty}\f$ and \f$f_{sub}^{uncertainty}\f$
 * are left to the uncertainty object to deal with, with the operators + and -.
 *
 * As uncertainty and unit are pointers and are replaceable, a memory management is
 * put into into place. A boolean tests if the program takes care of the deleting
 * of the pointer or if it is to the user to do so. By default, all is managed by
 * the program.
 * The object setters (void setUnitObject(Units*) and void setUncertaintyObject(CoreUnc*))
 * are the methods to use to switch to user management.
 *
 * The unit pointer has a special treatment, as densification can be perform. Thus,
 * the deleting process first supress one connection to the object, and delete it
 * only if there is no more ParameterPhy associated to it.
 *
 * Biaises in uncertainty calculations can happen: the program will calculate every uncertainty
 * in a serial way, which mean that any variable that appears several times
 * in a expression will be treated independently. For instance, consider the
 * expression \f$h(x) = f(x)g(x)\f$, the program will first compute the uncertainties
 * due to expression \f$g\f$ and combine it to the uncertainties due to expression
 * \f$f\f$, thus giving in the end:
 * \f[
 *    u_{h_c}^2 = g(x)^2  u_{f(x)}^2 + f(x)^2  u_{g(x)}^2
 * \f]
 * with \f$u_{f(x)}^2 = \left(\frac{\partial f}{\partial x}\right)^2  u_x^2\f$, and
 * the respective expression for \f$u_{g(x)}^2\f$.
 *
 * Yet the expression proposed by the GUM (and used in this program) yields:
 * \f[
 *      u_{h_{GU\!M}}^2 = \left[\left(\frac{\partial f}{\partial x}\right)  g(x) + 
 *                            \left(\frac{\partial g}{\partial x}\right)  f(x)  \right]^2  u_x^2
 * \f]
 * Those two expression differs by:
 * \f[
 *      biais_{u_h^2}^{mult} = u_{h_{GU\!M}}^2 - u_{h_c}^2 = 2  \left(\frac{\partial f}{\partial x}\right)  
 *                                                    \left(\frac{\partial g}{\partial x}\right) 
 *                                                    g(x)  f(x)  u_x^2
 * \f]
 * For a sum (\f$h(x) = f(x) + g(x)\f$), we obtain:
 * \f[
 *      biais_{u_h^2}^{add} = u_{h_{GU\!M}}^2 - u_{h_c}^2 = 2  \left(\frac{\partial f}{\partial x}\right)  
 *                                                    \left(\frac{\partial g}{\partial x}\right)  u_x^2
 * \f]
 */
class ParameterPhy{
  public:
/*!\brief Copy constructor, uses ParameterPhy &operator=(const ParameterPhy&)*/
   ParameterPhy(const ParameterPhy &ToBeCopied):
                unit(NULL),uncertainty(NULL),
                pdf(NULL),managePtrUnit(true),managePtrUnc(true)
                {*this = ToBeCopied;}
/*!\brief Simplest building constructor*/
   ParameterPhy(AtomicParameter par, Units *unitPtr, CoreUnc *uncPtr);
/*!\brief Constructor from AtomicParameter*/
   ParameterPhy(AtomicParameter par,std::string unitStr = std::string(),double dValue = -1.,std::string dUnitStr = CORE_UNCERTAINTY_TYPE_UNKNOWN);
/*!\brief Constructor from AtomicParameter*/
   ParameterPhy(AtomicParameter par,double dValue,std::string dUnitStr = CORE_UNCERTAINTY_TYPE_UNKNOWN, std::string unitStr = std::string());
/*!\brief Constructor from AtomicParameter*/
   ParameterPhy(AtomicParameter par,std::vector<double> dValues,std::string dUnitStr, std::string unitStr = std::string());
/*!\brief Constructor from AtomicParameter*/
   ParameterPhy(AtomicParameter par,std::string unitStr, std::vector<double> dValues,std::string dUnitStr);

/*!\brief User-friendly constructor
 *
 * Only a name and maybe a unit and an uncertainty type
 * are required.
 * Default values are no unit and unknown uncertainty
 */
   ParameterPhy(std::string namePar, std::string UnitStr = std::string(), std::string UncStr = CORE_UNCERTAINTY_TYPE_UNKNOWN);
/*!\brief Constructor overload, to cover possibilities*/
   ParameterPhy(std::string namePar, std::vector<double> valPar,
                std::string unitStr = std::string(),double dValue = -1., std::string dUnitStr = CORE_UNCERTAINTY_TYPE_UNKNOWN);
/*!\brief Constructor overload, to cover possibilities*/
   ParameterPhy(std::string namePar, std::vector<double> valPar,
                double dValue, std::string dUnitStr, std::string unitStr = std::string());
/*!\brief Constructor overload, to cover possibilities*/
   ParameterPhy(std::string namePar, std::vector<double> valPar,
                std::vector<double> dValues, std::string dUnitStr = CORE_UNCERTAINTY_TYPE_UNKNOWN, std::string unitStr = std::string());
/*!\brief Constructor overload, to cover possibilities*/
   ParameterPhy(std::string namePar, std::vector<double> valPar,
                std::string unitStr, std::vector<double> dValues, std::string dUnitStr);

/*!\brief Constructor to store only a value*/
  ParameterPhy(double value);

/*!\brief Constructor adapted to describe constant physical quantities*/
   ParameterPhy(std::string namePar,double valuePhy,double dvaluePhy,std::string typeDvaluePhy,std::string unitPar);

/*!\brief Default constructor*/
   ParameterPhy();
                  

/*!\brief Default destructor with memory management.
 *
 * If the memory management is automatic, then the deleting process
 * is launched:
 *    - the uncertainty memory is deallocated,
 *    - one connection is supressed from the unit object,
 *    if there is no more connection to it, the memory is
 *    deallocated.
 */
   ~ParameterPhy();

/*!\brief Name getter*/
   const std::string namePar()                   const {return parameter.getName();}
/*!\brief Name setter*/
         void  setNamePar(const std::string &parName)  {parameter.setName(parName);}

/********************************************
 *
 * Unit
 *
 * ******************************************/

/*!\brief Unit object getter*/
         Units *getUnitObject()                 const {return unit;}
/*!\brief Unit object setter*/
         void setUnitObject(Units *target);
/*!\brief Unit object changer
 *
 * This method is different from the setter in two ways:
 *    - it does not switch to user management,
 *    - the object given in the parameter is not considered
 *    to be initialized, it is initialized through the setUnit(std::string)
 *    method (which oblige the targetted unit object to have one).
 */
         void changeUnitObject(Units *target);
/*!\brief Create a copy of the given object, and assign it to the parameter*/
         void copyUnitObject(Units *target);
/*!\brief Unit densification
 *
 * Densification of unit means equalizing the pointers
 * if the units are homogeneous, with the appropriate
 * change in the values, if needed.
 * This method has the advantage to harmonize the
 * unit system used for all the ParameterPhy taken to be densified.
 *
 */
         bool densifyUnitObject(const ParameterPhy &targetPhy);
/*!\brief Advanced densification, the current object is changed to match the unit, if needed
 *
 * This is but an alias for bool densifyUnitObjectAlways(const ParameterPhy &),
 * here the idea is only to have a user-friendly name-explicit method
 *
 */
         bool densifyMeToThis(const ParameterPhy &targetPhy)    {return densifyUnitObject(targetPhy);}
/*!\brief Advanced densification, the given object is changed to match the unit, if needed
 *
 * This is the counter part of bool densifyMeToThis(const ParameterPhy &).
 */
         bool densifyThisToMe(ParameterPhy &targetPhy)    const {return targetPhy.densifyUnitObject(*this);}

/*!\brief Unit setter*/
         void setUnitPar(const std::string &parUnit);
/*!\brief Change the unit of the parameter to a wanted unit.
 *
 * If the units are not homogeneous, an error is returned, then,
 * if the unit object is a Units, nothing is done,
 * then the calculations are performed, on the uncertainty and
 * the values. It sets the uncertainty type to CORE_UNCERTAINTY_TYPE_ABSOLUTE.
 */
         void changeToSomeUnit(const std::string &target);
/*!\brief Change the unit of the parameter to the SI unit system.*/
         void changeToSIUnit()                           {changeToSomeUnit(unit->getSISymb());}
/*!\brief Change the unit of the parameter to the SI unit system, checking for SI composed unit.*/
         void changeToSIConvenientUnit()                 {changeToSomeUnit(unit->getSIConvenientSymb());}
/*!\brief Change the unit of the parameter to the harmonized symbol (see class Units for details).*/
         void harmonizeUnit();
/*!\brief Change the unit of the parameter to the contracted symbol (see class Units for details).*/
         void contractUnit();

/*!\brief Boolean to test the unit has a dimension.*/
         const bool isUnited()                               const;
/*!\brief Homogeneity test.*/
         const bool isHomogeneous(std::string target)        const;
/*!\brief Homogeneity test.*/
         const bool isHomogeneous(const ParameterPhy &rhs)   const {return isHomogeneous(rhs.unitPar());}

/*!\brief Unit getter.*/
         const std::string unitPar()                         const;
/*!\brief Returns the factor to a given unit. The unit given must be homogeneous to the current unit.*/
         double getFactorToSomeUnit(const std::string &target) const;
/*!\brief Returns the factor to the SI basis.*/
         double getFactor()                        const;
/*!\brief Returns the translator to a given unit. The unit given must be homogeneous to the current unit.*/
         double getTranslatorToSomeUnit(const std::string &target) const;
/*!\brief Returns the translator to the SI basis.*/
         double getTranslator()                    const;


/*****************************************
 *
 * Uncertainty
 * ***CoreUnc
 * ***************************************/


/*!\brief Uncertainty object getter*/
         CoreUnc *getUncertaintyObject()     const {return uncertainty;}
/*!\brief Uncertainty object setter*/
         void setUncertaintyObject(CoreUnc *target);

/*!\brief Uncertainty type getter*/
   const std::string typeDPar()                 const;
/*!\brief Uncertainty type setter*/
         void setTypeDPar(const std::string &dparType);

/*!\brief Uncertainty value getter.*/
         double Dvalue(int jval = 0)       const;
/*!\brief Overload of double Dvalue(int).*/
         double Dvalue(unsigned int jval)  const  {return Dvalue((int)jval);}
/*!\brief Variance value getter.*/
         double variance(int jval = 0)     const;
/*!\brief Overload of double variance(int).*/
         double variance(unsigned int jval)const  {return variance((int)jval);}
/*!\brief Standard uncertainty value getter.*/
         double StdUnc(int jval = 0)     const    {return (std::sqrt(variance(jval)));}
/*!\brief Overload of double StdUnc(int).*/
         double StdUnc(unsigned int jval)const    {return StdUnc((int)jval);}
/*!\brief Uncertainty value setter.*/
         void setDValue(double a,int i = 0);
/*!\brief Overload of double setDValue(int).*/
         void setDValue(double a,unsigned int i) {setDValue(a,(int) i);}
/*!\brief Variance value setter.*/
         void setVariance(double a,int i = 0);
/*!\brief Overload of double setVariance(int).*/
         void setVariance(double a,unsigned int i) {setVariance(a,(int) i);}
/*!\brief Standard uncertainty value setter.*/
         void setStdUnc(double a,int i = 0);
/*!\brief Overload of double setStdUnc(int).*/
         void setStdUnc(double a,unsigned int i) {setStdUnc(a,(int) i);}
/*!\brief Uncertainty value adder.*/
         void addDValue(double a);
/*!\brief Variance value adder.*/
         void addVariance(double var);
/*!\brief Standard uncertainty value adder.*/
         void addStdUnc(double stdU);
/*!\brief Correct an additive biais on uncertainty, on the variance*/
         void correcUncAddBiais(const ParameterPhy& correction,const ParameterPhy &delta);

/*!\brief Vector of uncertainty values getter.*/
   const std::vector<double> DValues()           const;
/*!\brief Vector of variances getter.
 *
 * If the uncertainty is of type CORE_UNCERTAINTY_TYPE_UNKNOWN
 * it sends back an empty vector, if the uncertainty is 
 * CORE_UNCERTAINTY_TYPE_NONE, it sends back only one value (0.)
 * in the vector. Else, it sends back the variance.
 */
   const std::vector<double> variances()         const;
/*!\brief Vector of uncertainty value setter.*/
         void setDValues(const std::vector<double> &a);
/*!\brief Vector of uncertainty value adder.*/
         void addDValues(const std::vector<double> &a);

/*!\brief Check if uncertainty is of type CORE_UNCERTAINTY_TYPE_ABSOLUTE*/
         bool absoluteUnc()                      const {return (uncertainty->getType() == CORE_UNCERTAINTY_TYPE_ABSOLUTE);}
/*!\brief Check if uncertainty is of type CORE_UNCERTAINTY_TYPE_RELATIVE*/
         bool relativeUnc()                      const {return (uncertainty->getType() == CORE_UNCERTAINTY_TYPE_RELATIVE);}
/*!\brief Check if uncertainty is of type CORE_UNCERTAINTY_TYPE_NONE*/
         bool noUnc()                            const {return (uncertainty->getType() == CORE_UNCERTAINTY_TYPE_NONE);}
/*!\brief Check if uncertainty is of type CORE_UNCERTAINTY_TYPE_UNKNOWN*/
         bool unknownUnc()                       const {return (uncertainty->getType() == CORE_UNCERTAINTY_TYPE_UNKNOWN);}
/*!\brief Check if uncertainty is other than CORE_UNCERTAINTY_TYPE_UNKNOWN or CORE_UNCERTAINTY_TYPE_NONE*/
         bool isUncDef()                         const {return (!unknownUnc() && !noUnc());}

/*!\brief Change the uncertainty type to CORE_UNCERTAINTY_TYPE_ABSOLUTE 
 *
 * Only performs the change from CORE_UNCERTAINTY_TYPE_RELATIVE
 * to CORE_UNCERTAINTY_TYPE_ABSOLUTE. Other cases need not to
 * be treated.
 */
         void setUncertaintyToAbsolute();
/*!\brief Change the uncertainty type to CORE_UNCERTAINTY_TYPE_RELATIVE
 *
 * Computes the necessary changes of values, does
 * nothing in case of CORE_UNCERTAINTY_TYPE_UNKNOWN
 */
         void setUncertaintyToRelative();
/*!\brief Provide the absolute uncertainty value, i.e. the variance.
 *
 *  If the type is UNCERTAINTY_TYPE_UNKNOWN or an unknown
 *  type, it returns an error. Else, the absolute value
 *  is computed and returned.
 *  */
         double getDValueAbsolute(int jval)const;
/*!\brief Overload of double getDValueAbsolute(int)const*/
         double getDValueAbsolute(unsigned int jval) const {return getDValueAbsolute((int)jval);}

/*!\brief Number of uncertainty values*/
         int nDValues()                     const;

/*****************************************
 *
 * Uncertainty
 * ***Pdf
 * ***************************************/

/*!\brief Return the name of the pdf (ex: Normal, Truncated normal,...)*/
   const std::string getPdfName()               const {return pdf->getNamePdf();}
/*!\brief Return the id of the pdf (ex: norm, nort,...)
 *
 * A pdf ID has two rules:
 *    - 4 characters
 *    - lower case
 * These two rules are to insure compatibility with a previous
 * fortran code from the coder, because I'm lazy and I don't
 * want to rewrite my input files.
 * */
   const std::string getPdfID()       const;
/*!\brief Returns the mean value for each pdf. See Pdf.cpp for detailed computation of each case.*/
         double getMean()        const;
/*!\brief Returns the median value for each pdf. See Pdf.cpp for detailed computation of each case.*/
         double getMedian()      const;
/*!\brief Returns the nominal value for each pdf. See Pdf.cpp for detailed computation of each case.*/
         double getMode()        const;
/*!\brief Returns the standard deviation for each pdf. See Pdf.cpp for detailed computation of each case.*/
         double getStandardDev()  const;
/*!\brief Returns the upper bound value if uniform, 4-sigma limit else*/
         double getUpperBound()  const;
/*!\brief Returns the lower bound value if uniform, 4-sigma limit else*/
         double getLowerBound()  const;


/*\brief Set the pointer of the pdf*/
         void setParameterPfd(Pdf *ptrPdf)  {pdf = ptrPdf;pdf->addOneConnection();}
/*\brief Set the pdf using the ID and the parameters
 *
 * No more than 4 par, -1e303 value chosen as never expected (\f$-\infty\f$), 
 * construction of std::vector<double> with the correct number of pars
 * so dealt by Pdf *objectPDF(const std::string &,const std::vector<double> &);
 *
 * */
         void setParameterPdf(const std::string &pdfID,
                                        double p1=-1e303,
                                        double p2=-1e303,
                                        double p3=-1e303,
                                        double p4=-1e303);
/*\brief Add this parameter to the pdf using the ID and the parameters
 *
 * No more than 4 par, -1e303 value chosen as never expected (\f$-\infty\f$), 
 * construction of std::vector<double> with the correct number of pars
 * */
         void addParameterToPdf(double p1=-1e303,
                                        double p2=-1e303,
                                        double p3=-1e303,
                                        double p4=-1e303);

         void addParameterToPdf(const std::vector<double> &data);

        void setParameterPdf(const std::string &pdfId,const std::vector<double> &data, int i = -1);
/*!\brief Test if the pdf is ok, this triggers the pdf tests*/
         void checkParameterPdf() {pdf->isItGood();}        


/*************************************************
 *
 * Values
 *
 * ***********************************************/


/*!\brief Number of values*/
         int nValues()                    const {return parameter.sizeVal();}
/*!\brief Value getter.*/
         double value(unsigned int jval)  const {return parameter.getValue(jval);}
/*!\brief Value getter.*/
         double value (int jval = 0)      const {return parameter.getValue(jval);}
/*!\brief First value getter*/
         double front()                   const {return parameter.front();}
/*!\brief Last value getter*/
         double back()                    const {return parameter.back();}
/*!\brief Value setter.*/
   void setValue(double a, int i = 0)            {parameter.setValue(a,i);}
/*!\brief Overload of void setValue(double,int).*/
   void setValue(double a, unsigned int i)       {setValue(a,(int)i);}
/*!\brief Value adder.*/
   void addValue(double a)                       {parameter.addValue(a);}

/*!\brief Vector of values setter*/
   void setValues(const std::vector<double> &a)        {parameter.setValues(a);}
/*!\brief Vector of values adder*/
   void addValues(const std::vector<double> &a)        {parameter.addValues(a);}
/*!\brief Vector of values getter*/
   const std::vector<double> values()            const {return parameter.getValues();}

/*\brief Add an already existing value a given number of times*/
   void duplicateValue(int nval, int ntimes)                   {parameter.duplicateVal(nval,ntimes);}
/*\brief Add an already existing value a given number of times, unsigned int overload*/
   void duplicateValue(unsigned int nval, unsigned int ntimes) {parameter.duplicateVal(nval,ntimes);}

/*****************************************
 *
 * from std::vector.h
 * 
 * ***************************************/


/*!\brief Resize the std::vector of values*/
   void resizePar(int newSize, double val = 0.)          {parameter.resizeValues(newSize,val);}
/*!\brief Resize the std::vector of values*/
   void resizePar(unsigned int newSize, double val = 0.) {parameter.resizeValues(newSize,val);}
/*!\brief Resize the std::vector of uncertainty values*/
   void resizeUnc(int newSize, double dv = 0.)           {uncertainty->resizeDValues(newSize,dv);}
/*!\brief Resize the std::vector of uncertainty values*/
   void resizeUnc(unsigned int newSize, double dv = 0.)  {uncertainty->resizeDValues(newSize,dv);}
/*!\brief Resize the std::vector of values and the uncertainty values*/
   void resize(int newSize, double val = 0., double dval = 0.)          {resizePar(newSize,val);resizeUnc(newSize,dval);}
/*!\brief Resize the std::vector of values and the uncertainty values*/
   void resize(unsigned int newSize, double val = 0., double dval = 0.) {resizePar(newSize,val);resizeUnc(newSize,dval);}
/*!\brief Reserve the std::vector of values*/
   void reservePar(int newSize)                       {parameter.reserveValues(newSize);}
/*!\brief Reserve the std::vector of values*/
   void reservePar(unsigned int newSize)              {parameter.reserveValues(newSize);}
/*!\brief Reserve the std::vector of uncertainty values*/
   void reserveUnc(int newSize)                       {uncertainty->reserveDValues(newSize);}
/*!\brief Reserve the std::vector of uncertainty values*/
   void reserveUnc(unsigned int newSize)              {uncertainty->reserveDValues(newSize);}
/*!\brief Reserve the std::vector of values and the uncertainty values*/
   void reserve(int newSize)                         {reservePar(newSize);reserveUnc(newSize);}
/*!\brief Reserve the std::vector of values and the uncertainty values*/
   void reserve(unsigned int newSize)                {reservePar(newSize);reserveUnc(newSize);}


/**************************************************
 *
 * Global
 *
 * ************************************************/

/*!\brief Replace the parameter
 *
 * This is the full copy, all is equalized with
 * the ParameterPhy &operator=(const ParameterPhy &), and
 * then the name is changed. This is not equality, this is
 * replacement.
 */
   void replace(const ParameterPhy &rhs);
/*!\brief ParameterPhy adder.
 *
 * Adding the values and uncertainty values of
 * another parameter to the current values and
 * uncertainty values needs several checking and
 * some processing.
 *
 * A decision is which parameter ajusts to the other.
 *
 * Finally, the uncertainty are set to absolute and
 * the values are added.
 *
 * This version is the default version: advanced unit
 * management and adjustment of the given parameter.
 */
   void complete(const ParameterPhy &comp);

/*!\brief AtomicParameter getter*/
   const AtomicParameter getParameter()      const {return parameter;}

/*!\brief Is equal operator. Two physical parameters are equal if their names are the same.*/
   bool const     operator==(const ParameterPhy &rhs) const {return (namePar() == rhs.namePar());}
/*!\brief Not is equal operator.*/
   bool const     operator!=(const ParameterPhy &rhs) const {return (!(*this == rhs));}
/*!\brief Assignement operator. Uses the assignement operator of units, uncertainty and atomic parameter.*/
   ParameterPhy & operator= (const ParameterPhy &rhs);
/*!\brief Adding operator. 
 *
 * Uses the adding operator of uncertainty and atomic parameter.
 * If the units are not equal, the operator calls the void changeToSomeUnit(std::string)
 * method.
 */
   ParameterPhy & operator+=(const ParameterPhy &rhs);
/*!\brief Substracting operator. This is adding the parameter multiplied by the double -1.*/
   ParameterPhy & operator-=(const ParameterPhy &rhs);
/*!\brief Multiplicating operator with a double.
 *
 * Changes the values and uncertainty values if the
 * uncertainty type is CORE_UNCERTAINTY_TYPE_ABSOLUTE
 */
   ParameterPhy & operator*=(double r);
/*!\brief Dividing operator with a double. See ParameterPhy & operator*=(double)*/
   ParameterPhy & operator/=(double r);
/*!\brief Adding operator with a double.
 *
 * Adding a scalar is a translation, so no unit consideration is needed,
 * The absolute uncertainty does not change. To prevent biais in case
 * of relative uncertainty, it is set to absolute.
 */
         ParameterPhy & operator+=(double r);
/*!\brief Substracting operator with a double. This is adding the opposite.*/
         ParameterPhy & operator-=(double r);
/*!\brief Multiplicating operator. 
 *
 * Uses the adding operator of the unit, the multiplying operator of the
 * AtomicParameter object, and the GUM expression for the uncertainty:
 * \f[
 *      f_{mult}^{uncertainty}(\varphi_1 \times \varphi_2) = u_c^2 = \nu_2^2u_1^2 + \nu_1^2u_2^2
 * \f]
 * with \f$\nu_i\f$ being the value of parameter \f$i\f$ and \f$u\f$ the standard
 * uncertainty.
 */
         ParameterPhy & operator*=(const ParameterPhy &rhs);
/*!\brief Dividing operator.
 *
 * See ParameterPhy & operator*=(const ParameterPhy &);
 * The uncertainty model for division is:
 * \f[
 *     f_{div}^{uncertainty}(\frac{\varphi_1}{\varphi_2}) = 
 *     u_c^2 = \left(\frac{1}{\nu_2}\right)^2u_1^2 + \left(-\frac{\nu_1}{\nu_2^2}\right)^2u_2^2
 * \f]
 */
         ParameterPhy & operator/=(const ParameterPhy &rhs);
/*!\brief Multiplicating operator. Uses the multiplicating with a double operator.*/
         ParameterPhy   operator* (double r) const;
/*!\brief Dividing operator.  Uses the dividing with a double operator.*/
         ParameterPhy   operator/ (double r) const;
/*!\brief Adding operator.  Uses the additing with a double operator.*/
         ParameterPhy   operator+ (double r) const;
/*!\brief Substracting operator.  Uses the substracting with a double operator.*/
         ParameterPhy   operator- (double r) const;
/*!\brief Multiplicating operator.  Uses the multiplicating with a double operator.*/
         ParameterPhy   operator* (const ParameterPhy &rhs) const;
/*!\brief Dividing operator.  Uses the dividing with a double operator.*/
         ParameterPhy   operator/ (const ParameterPhy &rhs) const;
/*!\brief Adding operator.  Uses the adding with a double operator.*/
         ParameterPhy   operator+ (const ParameterPhy &rhs) const;
/*!\brief Substracting operator.  Uses the substracting.*/
         ParameterPhy   operator- (const ParameterPhy &rhs) const;
/*!\brief Inversing the values, nothing else changes*/
         ParameterPhy   operator- () const {return (ParameterPhy(-parameter,unit,uncertainty));}
/*!\brief Friend operator to automatize calculations with doubles*/
         friend ParameterPhy operator+ (const double &lhs, const ParameterPhy &rhs);
/*!\brief Friend operator to automatize calculations with doubles*/
         friend ParameterPhy operator- (const double &lhs, const ParameterPhy &rhs);
/*!\brief Friend operator to automatize calculations with doubles*/
         friend ParameterPhy operator* (const double &lhs, const ParameterPhy &rhs);
/*!\brief Friend operator to automatize calculations with doubles*/
         friend ParameterPhy operator/ (const double &lhs, const ParameterPhy &rhs);

/*!\brief Alternative to ParameterPhy &operator=(const ParameterPhy&)*/
         void equalize(const ParameterPhy &target)  {*this = target;}
/*!\brief Alternative to ParameterPhy &operator+=(const ParameterPhy&)*/
         void add(const ParameterPhy &target)       {*this += target;}
/*!\brief Alternative to ParameterPhy &operator+=(double)*/
         void add(double r)                         {*this += r;}
/*!\brief Alternative to ParameterPhy &operator-=(const ParameterPhy&)*/
         void substract(const ParameterPhy &target) {*this -= target;}
/*!\brief Alternative to ParameterPhy &operator-=(double)*/
         void substract(double r)                   {*this -= r;}
/*!\brief Alternative to ParameterPhy &operator*=(const ParameterPhy&)*/
         void multiply(const ParameterPhy &target)  {*this *= target;}
/*!\brief Alternative to ParameterPhy &operator*=(double)*/
         void multiply(double r)                    {*this *= r;}
/*!\brief Alternative to ParameterPhy &operator/=(const ParameterPhy&)*/
         void divide(const ParameterPhy &target)    {*this /= target;}
/*!\brief Alternative to ParameterPhy &operator/=(double)*/
         void divide(double r)                      {*this /= r;}


/*!\brief Clear the parameter values*/
  void clear();
/*!\brief Test if the parameter is empty or not*/
  bool empty()      const {return (parameter.sizeVal() == 0);}
/*!\brief Clear the parameter's values and name*/
  void clearParameter()   {parameter.clear();}
/*!\brief Clear the parameter's values*/
  void clearValues()      {parameter.clearValues();}
/*!\brief Clear the unit values*/
  void clearUnit()        {unit->clear();}
/*!\brief Clear the uncertainty values*/
  void clearUncertainty() {uncertainty->clear();}

/*!\brief Expose the parameter*/
  void showAll(std::ostream &out = std::cout) const;

/*!\brief managePtrUnit getter*/
  bool getManagementUnit(){return managePtrUnit;}
/*!\brief managePtrUnc getter*/
  bool getManagementUnc() {return managePtrUnc;}

/*!\brief managePtrUnit setter to manual unit managment*/
  void IManageTheUnit(){managePtrUnit = false;}
/*!\brief managePtrUnc setter to manual uncertainty managment*/
  void IManagetheUnc() {managePtrUnc = false;}

  private:

/*!\brief ParameterPhy adder root for the version that uses a const ParameterPhy& as parameter.*/
   void completeRoot(const ParameterPhy &comp,bool advanced);
/*!\brief GUM advised model for the combined uncertainty, applied to multiplication.
 *
 * \f[
 *      f_{add}^{uncertainty}(\varphi_1 + \varphi_2) = u_c^2 = u_1^2 + u_2^2
 * \f]
 */
  double additiveModelUncertainty(const double &val1, const double &val2, const double &var1, const double &var2);
/*!\brief GUM advised model for the combined uncertainty, applied to multiplication.
 *
 * \f[
 *      f_{mult}^{uncertainty}(\varphi_1 \times \varphi_2) = u_c^2 = \nu_2^2u_1^2 + \nu_1^2u_2^2
 * \f]
 */
  double multiplyModelUncertainty(const double &val1, const double &val2, const double &var1, const double &var2);
/*!\brief GUM advised model for the combined uncertainty, applied to division.
 *
 * \f[
 *     f_{div}^{uncertainty}(\frac{\varphi_1}{\varphi_2}) = 
 *     u_c^2 = \left(\frac{1}{\nu_2}\right)^2u_1^2 + \left(-\frac{\nu_1}{\nu_2^2}\right)^2u_2^2
 * \f]
 */
  double divideModelUncertainty(const double &val1, const double &val2, const double &var1, const double &var2);
/*!\brief Uncertainty manager
 *
 * It will combine the uncertainty of two parameters, using the given
 * model, set to the function pointer modelU(double value1, double value2, double variance1, double variance2).
 * double divideModelUncertainty(double val1, double val2, double var1, double var2) and
 * double multiplyModelUncertainty(double val1, double val2, double var1, double var2)
 * are exemples.
 */
  void modelUncertainty(const ParameterPhy &rhs,double (ParameterPhy::*modelU)(const double&,const double&,const double&,const double&));

/*! \brief Set the pointer to NULL value, with memory management.
 *
 * The memory is deallocated only if this parameter is the only
 * one connected to the unit object pointed by the pointer.
 */
  void UnitPtrToNull();

/*!\brief An AtomicParameter for the parameter.*/
   AtomicParameter parameter;
/*!\brief A Units pointer for the unit.*/
   Units *unit;
/*!\brief A CoreUnc pointer for the uncertainty.*/
   CoreUnc *uncertainty;
/*!\brief A Pdf pointer for the pdf*/
   Pdf *pdf;
/*!\brief Two booleans for the memory management of the unit and uncertainty pointers.*/
   bool managePtrUnit,managePtrUnc;

};

}
#endif
