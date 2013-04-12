//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

//Antioch
#include "antioch/ParameterPhy.hpp"

//C++

namespace Antioch{

ParameterPhy::ParameterPhy():
        unit(NULL),uncertainty(NULL),
        pdf(NULL),managePtrUnit(true),
        managePtrUnc(true)
{
  unit = new Units;
  unit->addOneToConnection();
  uncertainty = new CoreUnc(CORE_UNCERTAINTY_TYPE_UNKNOWN);
}


ParameterPhy::ParameterPhy(AtomicParameter par, Units *unitPtr, CoreUnc *uncPtr):
              unit(NULL),uncertainty(NULL),pdf(NULL),
              managePtrUnit(true),managePtrUnc(true)
{
  parameter = par;
  copyUnitObject(unitPtr);
  uncertainty = new CoreUnc(*uncPtr);
}

ParameterPhy::ParameterPhy(AtomicParameter par,std::string unitStr,double dValue,std::string dUnitStr):
              parameter(par),
              unit(NULL),uncertainty(NULL),pdf(NULL),
              managePtrUnit(true),managePtrUnc(true)
{
  unit = new Units(unitStr);
  unit->addOneToConnection();
  uncertainty = new CoreUnc(dValue,dUnitStr);
}

ParameterPhy::ParameterPhy(AtomicParameter par,double dValue,std::string dUnitStr, std::string unitStr):
              parameter(par),
              unit(NULL),uncertainty(NULL),pdf(NULL),
              managePtrUnit(true),managePtrUnc(true)
{
  unit = new Units(unitStr);
  unit->addOneToConnection();
  uncertainty = new CoreUnc(dValue,dUnitStr);
}

ParameterPhy::ParameterPhy(AtomicParameter par,std::vector<double> dValues,std::string dUnitStr, std::string unitStr):
              parameter(par),
              unit(NULL),uncertainty(NULL),pdf(NULL),
              managePtrUnit(true),managePtrUnc(true)
{
  unit = new Units(unitStr);
  unit->addOneToConnection();
  uncertainty = new CoreUnc(dValues,dUnitStr);
}

ParameterPhy::ParameterPhy(AtomicParameter par,std::string unitStr, std::vector<double> dValues,std::string dUnitStr):
              parameter(par),
              unit(NULL),uncertainty(NULL),pdf(NULL),
              managePtrUnit(true),managePtrUnc(true)
{
  unit = new Units(unitStr);
  unit->addOneToConnection();
  uncertainty = new CoreUnc(dValues,dUnitStr);
}

//use-friendly, same thing
ParameterPhy::ParameterPhy(std::string namePar, std::string UnitStr, std::string UncStr):
              parameter(namePar),
              unit(NULL),uncertainty(NULL),pdf(NULL),
              managePtrUnit(true),managePtrUnc(true)
{
  unit = new Units(UnitStr);
  unit->addOneToConnection();
  if(UncStr.empty())UncStr = CORE_UNCERTAINTY_TYPE_UNKNOWN;
  uncertainty = new CoreUnc(UncStr);
}


ParameterPhy::ParameterPhy(std::string namePar, std::vector<double> valPar, std::string unitStr,double dValue, std::string dUnitStr):
              parameter(namePar,valPar),
              unit(NULL),uncertainty(NULL),pdf(NULL),
              managePtrUnit(true),managePtrUnc(true)
{
  unit = new Units(unitStr);
  unit->addOneToConnection();
  uncertainty = new CoreUnc(dValue,dUnitStr);
}

ParameterPhy::ParameterPhy(std::string namePar, std::vector<double> valPar,double dValue, std::string dUnitStr, std::string unitStr):
              parameter(namePar,valPar),
              unit(NULL),uncertainty(NULL),pdf(NULL),
              managePtrUnit(true),managePtrUnc(true)
{
  unit = new Units(unitStr);
  unit->addOneToConnection();
  uncertainty = new CoreUnc(dValue,dUnitStr);
}

ParameterPhy::ParameterPhy(std::string namePar, std::vector<double> valPar,std::vector<double> dValues, std::string dUnitStr, std::string unitStr):
              parameter(namePar,valPar),
              unit(NULL),uncertainty(NULL),pdf(NULL),
              managePtrUnit(true),managePtrUnc(true)
{
  unit = new Units(unitStr);
  unit->addOneToConnection();
  uncertainty = new CoreUnc(dValues,dUnitStr);
}

ParameterPhy::ParameterPhy(std::string namePar, std::vector<double> valPar,std::string unitStr, std::vector<double> dValues, std::string dUnitStr):
              parameter(namePar,valPar),
              unit(NULL),uncertainty(NULL),pdf(NULL),
              managePtrUnit(true),managePtrUnc(true)
{
  unit = new Units(unitStr);
  unit->addOneToConnection();
  uncertainty = new CoreUnc(dValues,dUnitStr);
}

//constructor for a value
ParameterPhy::ParameterPhy(double valuePhy):
              parameter("scalar",valuePhy),
              unit(NULL),uncertainty(NULL),pdf(NULL),
              managePtrUnit(true),managePtrUnc(true)
{
  unit = new Units("");
  unit->addOneToConnection();
  uncertainty = new CoreUnc(CORE_UNCERTAINTY_TYPE_NONE);
}

ParameterPhy::~ParameterPhy()
{
  if(unit != NULL)
  {
    unit->suppressOneToConnection();
    if(managePtrUnit && unit->nConnections() == 0)delete unit;
  }
  if(uncertainty != NULL && managePtrUnc)delete uncertainty;
  if(pdf != NULL)
  {
    pdf->supressOneConnection();
    if(pdf->nConnection() <= 0)delete pdf;
  }
}

//constructor adapted for physical constant
ParameterPhy::ParameterPhy(std::string namePar,double valuePhy,double dvaluePhy,std::string typeDvaluePhy,std::string unitPar):
              parameter(namePar,valuePhy),
              unit(NULL),uncertainty(NULL),pdf(NULL),
              managePtrUnit(true),managePtrUnc(true)
{
  unit = new Units(unitPar);
  unit->addOneToConnection();
  uncertainty = new CoreUnc(dvaluePhy,typeDvaluePhy);
}

const std::vector<double> ParameterPhy::DValues() const
{
  return uncertainty->getDvalues();
}


//unit object
void ParameterPhy::changeUnitObject(Units *target)
{
  if(unit == target)return;
  target->setUnit(unit->getSymbol());
  UnitPtrToNull();
  unit = target;
  unit->addOneToConnection();
}

void ParameterPhy::setUnitObject(Units *target)
{
  UnitPtrToNull();
  if(target != NULL)
  {
    unit = target;
    unit->addOneToConnection();
  }
  managePtrUnit = false;
}

void ParameterPhy::copyUnitObject(Units *target)
{ 
  UnitPtrToNull();
  if(target == NULL)return;
  unit = new Units(target->getSymbol());
  unit->addOneToConnection();
  managePtrUnit = true;
}

bool ParameterPhy::densifyUnitObject(const ParameterPhy &targetPhy)
{
  Units *target = targetPhy.getUnitObject();
  if(target == NULL)return false;//if NULL ptr
  if(unit == target)return true;
  if(unit->isHomogeneous(target->getSymbol()))
  {
    changeToSomeUnit(target->getSymbol());
    setUnitObject(target);
    managePtrUnit = true;
    return true;
  }
  return false;
}

const bool ParameterPhy::isUnited() const
{
  return unit->isUnited();
}

const bool ParameterPhy::isHomogeneous(std::string target) const
{
  return unit->isHomogeneous(target);
}

const std::string ParameterPhy::unitPar() const
{ 
  return unit->getSymbol();
}

void ParameterPhy::setUnitPar(const std::string &parUnit)
{
  unit->setUnit(parUnit);
}

double ParameterPhy::getFactorToSomeUnit(const std::string &target) const
{
  return unit->FactorToSomeUnit(target);
}

double ParameterPhy::getTranslatorToSomeUnit(const std::string &target) const
{
  return unit->TranslatorToSomeUnit(target);
}

double ParameterPhy::getFactor() const
{
  return unit->getSIFactor();
}

double ParameterPhy::getTranslator() const
{
  return unit->getSITranslator();
}

void ParameterPhy::changeToSomeUnit(const std::string &target)
{
  if(unitPar() == target)return;
  if(!(unit->isHomogeneous(target)))
  {
    const std::string method("ParameterPhy::changeToSomeUnit(const std::string &)");
    std::string errStr = "\"" + unit->getSymbol() + 
                    "\" and \"" + target + "\" are not homogeneous!!\nNothing to do";
    antiochError(method,errStr);
  }else
  {
    double t = getTranslatorToSomeUnit(target);
    double v = getFactorToSomeUnit(target);
    if(isUncDef())
    {
      setUncertaintyToAbsolute();
      for(int i = 0; i < nDValues(); i++)
      {
        setDValue(v*Dvalue(i),i);
      }
    }
    for(int i = 0; i < nValues(); i++)
    {
      setValue(v*value(i) + t,i);
    }
    std::string tmp(target);
    unit->contractedSymbol(tmp);
    unit->setUnit(tmp);
  }
}

void ParameterPhy::contractUnit()
{
  setUnitPar(unit->contractedSymbol());
}

void ParameterPhy::harmonizeUnit()
{
  changeToSomeUnit(unit->harmonizedSymbol());
}

//uncertainty object
void ParameterPhy::setUncertaintyObject(CoreUnc *target)
{
  delete uncertainty;
  uncertainty = target;
}

//uncertainty type
const std::string ParameterPhy::typeDPar() const
{
  return uncertainty->getType();
}

void ParameterPhy::setTypeDPar(const std::string &dparType)
{
  uncertainty->setType(dparType);
}


//dvalue
double ParameterPhy::Dvalue(int jval) const
{
 return uncertainty->getDvalue(jval);
}

double ParameterPhy::variance(int jval) const
{
  if(typeDPar() == CORE_UNCERTAINTY_TYPE_ABSOLUTE || typeDPar() == CORE_UNCERTAINTY_TYPE_NONE)
  {
    return (uncertainty->getDvalue(jval));
  }else if(typeDPar() == CORE_UNCERTAINTY_TYPE_RELATIVE)
  {
    return (uncertainty->getDvalue(jval)*uncertainty->getDvalue(jval)*value(jval)*value(jval));
  }else
  {
    const std::string method("double ParameterPhy::variance(int) const");
    std::string err = "Cannot produce a variance with uncertainty type " + typeDPar();
    antiochError(method,err);
    return 0.;
  }
}

void ParameterPhy::setDValue(double a,int i)
{
  uncertainty->setDvalue(a,i);
}

void ParameterPhy::setVariance(double var,int jval)
{
  if(typeDPar() == CORE_UNCERTAINTY_TYPE_ABSOLUTE)
  {
    uncertainty->setDvalue(var,jval);
  }else if(typeDPar() == CORE_UNCERTAINTY_TYPE_RELATIVE)
  {
    uncertainty->setDvalue(sqrt(var)/value(jval),jval);
  }else
  {
    const std::string method("void ParameterPhy::setVariance(double,int) const");
    std::string err = "Cannot set a variance with uncertainty type " + typeDPar();
    antiochError(method,err);
  }
}

void ParameterPhy::setStdUnc(double stdU,int jval)
{
  if(typeDPar() == CORE_UNCERTAINTY_TYPE_ABSOLUTE)
  {
    uncertainty->setDvalue(stdU*stdU,jval);
  }else if(typeDPar() == CORE_UNCERTAINTY_TYPE_RELATIVE)
  {
    uncertainty->setDvalue(stdU/value(jval),jval);
  }else
  {
    const std::string method("void ParameterPhy::setStdUnc(double,int) const");
    antiochError(method,"Cannot set a standard uncertainty with uncertainty type " + typeDPar());
  }
}

void ParameterPhy::addDValue(double a)
{
  uncertainty->addDvalue(a);
}

void ParameterPhy::addVariance(double var)
{
  if(typeDPar() == CORE_UNCERTAINTY_TYPE_ABSOLUTE)
  {
    uncertainty->addDvalue(var);
  }else if(typeDPar() == CORE_UNCERTAINTY_TYPE_RELATIVE)
  {
    uncertainty->addDvalue(sqrt(var)/value(nDValues()));
  }else
  {
    const std::string method("void ParameterPhy::addVariance(double)");
    std::string err = "Cannot add a variance with uncertainty type " + typeDPar();
    antiochError(method,err);
  }
}

void ParameterPhy::addStdUnc(double stdU)
{
  if(typeDPar() == CORE_UNCERTAINTY_TYPE_ABSOLUTE)
  {
    uncertainty->addDvalue(sqrt(stdU));
  }else if(typeDPar() == CORE_UNCERTAINTY_TYPE_RELATIVE)
  {
    uncertainty->addDvalue(stdU/value(nDValues()));
  }else
  {
    const std::string method("void ParameterPhy::addStdUnc(double)");
    antiochError(method,"Cannot add a standard uncertainty with uncertainty type " + typeDPar());
  }
}

//dvalues
void ParameterPhy::setDValues(const std::vector<double> &a)
{
  uncertainty->setDvalues(a);
}

//dvalues
void ParameterPhy::addDValues(const std::vector<double> &a)
{
  uncertainty->addDvalues(a);
}

const std::string ParameterPhy::getPdfID() const
{
  const std::string method("const std::string ParameterPhy::getPdfID() const");

  if(pdf == NULL)antiochError(method,"Cannot produce an ID from a non-existant pdf.");

  return pdf->getIdPdf();
}

double ParameterPhy::getMean()        const
{
  const std::string method("double ParameterPhy::getMean() const");

  if(pdf == NULL)antiochError(method,"Cannot produce a mean from a non-existant pdf.");

  return pdf->getMeanMar();
}

double ParameterPhy::getMedian()      const
{
  const std::string method("double ParameterPhy::getMedian() const");

  if(pdf == NULL)antiochError(method,"Cannot produce a median from a non-existant pdf.");

  return pdf->getMedianMar();
}

double ParameterPhy::getMode()        const
{
  const std::string method("double ParameterPhy::getMode() const");

  if(pdf == NULL)antiochError(method,"Cannot produce a mode from a non-existant pdf.");

  return pdf->getModeMar();
}

double ParameterPhy::getStandardDev() const
{
  const std::string method("double ParameterPhy::getStandardDev() const");

  if(pdf == NULL)antiochError(method,"Cannot produce a standard deviation from a non-existant pdf.");

  return pdf->getStdDevMar();
}

double ParameterPhy::getUpperBound()  const
{
  const std::string method("double ParameterPhy::getUpperBound() const");

  if(pdf == NULL)antiochError(method,"Cannot produce an upper bound from a non-existant pdf.");

  return pdf->getUpLimitMar();
}

double ParameterPhy::getLowerBound()  const
{
  const std::string method("double ParameterPhy::getLowerBound() const");

  if(pdf == NULL)antiochError(method,"Cannot produce a lower bound from a non-existant pdf.");

  return pdf->getDownLimitMar();
}


void ParameterPhy::addParameterToPdf(double p1,double p2,double p3,double p4)
{
  std::vector<double> pars;
  if(p1 != -1e303)pars.push_back(p1);
  if(p2 != -1e303)pars.push_back(p2);
  if(p3 != -1e303)pars.push_back(p3);
  if(p4 != -1e303)pars.push_back(p4);
  
  addParameterToPdf(pars); 
}


void ParameterPhy::addParameterToPdf(const std::vector<double>&data)
{
  const std::string method("void ParameterPhy::addParameterToPdf(const std::vector<double>&)");
  if(pdf == NULL)
    antiochError(method,"Cannot add a parameter to a non existing pdf.");
  pdf->addMarginalFromData(data);
  pdf->addOneConnection();
}

void ParameterPhy::setParameterPdf(const std::string &pdfID,double p1,double p2,double p3,double p4)
{
  std::vector<double> pars;
  if(p1 != -1e303)pars.push_back(p1);
  if(p2 != -1e303)pars.push_back(p2);
  if(p3 != -1e303)pars.push_back(p3);
  if(p4 != -1e303)pars.push_back(p4);

  if(pdf != NULL)delete pdf;
  setParameterPdf(pdfID,pars);
}

void ParameterPhy::setParameterPdf(const std::string &pdfID,const std::vector<double> &data, int i)
{
  pdf = objectPDF(pdfID,data);
  if(i == -1)i=0;
  pdf->setNameMar(namePar(),i);
}

void ParameterPhy::modelUncertainty(const ParameterPhy &rhs,double (ParameterPhy::*modelU)(double,double,double,double))
{
  const std::string method("void ParameterPhy::modelUncertainty(const ParameterPhy &,double (ParameterPhy::*)(double, double, double, double))");

// if rhs has no unc object
  if(rhs.getUncertaintyObject() == NULL)
  {
     uncertainty->setTypeToUnknown();
     return;
  }

// unknown purged and contagious
  if(rhs.unknownUnc() || unknownUnc())
  {
     uncertainty->setTypeToUnknown();
     return;
  }

// now we go
  setUncertaintyToAbsolute();
  CoreUnc *mul = rhs.getUncertaintyObject();
  bool delmul(false);
  if(mul->getType() == CORE_UNCERTAINTY_TYPE_RELATIVE)
  {
     mul = NULL;
     delmul = true;
     mul = new CoreUnc(rhs.typeDPar());
     mul->setType(CORE_UNCERTAINTY_TYPE_ABSOLUTE);
     for(int i = 0; i < rhs.nDValues(); i++)
     {
        mul->addDvalue(rhs.getDValueAbsolute(i));
     }
  }
//now getting number of values
// hypothesis, nDValues() = nValues(), if not, there's a problem
  int min = (nValues() < rhs.nValues())?nValues():rhs.nValues();
  int max = (nValues() > rhs.nValues())?nValues():rhs.nValues();
  if((min != 1) && (max != min))
        antiochError(method,"No fixed parameter found (1 value), but the number of values are different for parameters "
                       + namePar() + " and " + rhs.namePar());
 
//setting the indexes, either from 0 to max - 1 or 0
  int i, z(0);
  int *kv1(&i),*kv2(&i),*ku1(&i),*ku2(&i);
  if(nValues() != max)kv1 = &z;
  if(nDValues() != max)
  {
    resizeUnc(max,Dvalue(0));
  }
  if(rhs.nValues() != max)kv2 = &z;
  if(rhs.nDValues() != max)ku2 = &z;
  for(i = 0; i < max; i++)
  {
    setDValue((this->*modelU)(value(*kv1),rhs.value(*kv2),Dvalue(*ku1),mul->getDvalue(*ku2)),i);
  }
  if(delmul)delete mul;
}

double ParameterPhy::additiveModelUncertainty(double val1, double val2, double var1, double var2)
{
  return (var2 + var1);
}

double ParameterPhy::multiplyModelUncertainty(double val1, double val2, double var1, double var2)
{
  return (val1*val1*var2 + val2*val2*var1);
}

double ParameterPhy::divideModelUncertainty(double val1, double val2, double var1, double var2)
{
  return (var1/(val2*val2) + var2*val1*val1/(val2*val2*val2*val2));
}

// conversion of dvalues in "absolute" core version
void ParameterPhy::setUncertaintyToAbsolute()
{
  if(typeDPar() != CORE_UNCERTAINTY_TYPE_ABSOLUTE)return;
  if(typeDPar() == CORE_UNCERTAINTY_TYPE_RELATIVE) // compute
  {
    double dv = uncertainty->getDvalue(0);
    bool update(nDValues() == nValues());
    for(int i = 0; i < nValues(); i++)
    {
      if(update)dv=uncertainty->getDvalue(i);
      setDValue(dv*value(i)*dv*value(i),i); 
    }
  }
  setTypeDPar(CORE_UNCERTAINTY_TYPE_ABSOLUTE);
  
}

// conversion of dvalues in "relative" core version
void ParameterPhy::setUncertaintyToRelative()
{
  if(typeDPar() == CORE_UNCERTAINTY_TYPE_ABSOLUTE) // compute
  {
    setTypeDPar(CORE_UNCERTAINTY_TYPE_RELATIVE);// relative protocol for getters => adjust
    for(int i = 0; i < nDValues(); i++)
    {
      setDValue(fabs(sqrt(Dvalue(i))/value(i)) , i); 
    }
  }else if(typeDPar() == CORE_UNCERTAINTY_TYPE_NONE)
  {
    uncertainty->clear();
    setTypeDPar(CORE_UNCERTAINTY_TYPE_RELATIVE);
    resizeUnc(nValues(),0.);
  }

}

double ParameterPhy::getDValueAbsolute(int jval) const 
{
  const std::string method("double ParameterPhy::getDValueAbsolute(int) const");
  if(typeDPar() == CORE_UNCERTAINTY_TYPE_ABSOLUTE) // already
  {
    return Dvalue(jval);
  }
  else if(typeDPar() == CORE_UNCERTAINTY_TYPE_RELATIVE) // compute
  {
    return fabs(Dvalue(jval)*value(jval)*Dvalue(jval)*value(jval));
  }
  else if(typeDPar() == CORE_UNCERTAINTY_TYPE_UNKNOWN) // no defaults, negative value returned
  {
    std::string err = "Parameter " + namePar() + ": cannot provide absolute uncertainty if type is CORE_UNCERTAINTY_TYPE_UNKNOWN.";
    antiochError(method,err);
    return -1.;
  }
  else if(typeDPar() == CORE_UNCERTAINTY_TYPE_NONE) // 0.
  {
    return 0.;
  }
  else // error
  {
    std::string err = "Parameter " + namePar() + ": cannot provide absolute uncertainty, type is unknown.";
    antiochError(method,err);
    return -1.;
  }
}

void ParameterPhy::correcUncAddBiais(const ParameterPhy& correction,const ParameterPhy &delta)
{
  const std::string method("void ParameterPhy::correcUncAddBiais(const ParameterPhy&, const ParameterPhy&)");
  if(correction.nValues() != nDValues() || delta.nDValues() != nDValues())
        {antiochError(method,"The correction and the uncertainty do not have the same number of values.");}

  setUncertaintyToAbsolute();
  for(int i = 0; i < correction.nValues(); i++)
  {
     setDValue(Dvalue(i) + correction.value(i)*delta.variance(i),i);
  }
}

//number of values, val & dval
int ParameterPhy::nDValues() const
{
  return uncertainty->sizeDVal();
}

void ParameterPhy::UnitPtrToNull()
{
  if(unit == NULL)return;
  unit->suppressOneToConnection();
  if(unit->nConnections() == 0)delete unit;
  unit = NULL;
}


void ParameterPhy::clear()
{
  parameter.clear();
  uncertainty->clear();
  unit->clear();
}

/*!\todo add pdf*/
void ParameterPhy::showAll(std::ostream &out) const
{
  out << namePar() << " in " << unitPar() << " with " << typeDPar() << " uncertainty" << std::endl;
  double dv = 0.;
  out << "  ";
  for(int i = 0; i < nValues(); i++)
  {
    if(isUncDef())dv = Dvalue(i);
    out << value(i) << " (" << dv << ") ";
  }
  out << std::endl;
}

void ParameterPhy::complete(const ParameterPhy &comp)
{
// test unit
  if(!isHomogeneous(comp))
  {
    const std::string method("void ParameterPhy::complete(const ParameterPhy &)");
    antiochError(method,"Parameters " + namePar() + " and " + comp.namePar() + " are not homogeneous");
  }

  double fact = comp.getFactorToSomeUnit(unitPar());
  if(!comp.unknownUnc() && !unknownUnc())
  {
    if(noUnc())
    {
       uncertainty->clear();
       uncertainty->resizeDValues(nValues(),0.);
       uncertainty->setType(CORE_UNCERTAINTY_TYPE_ABSOLUTE);
    }
    for(int i = 0; i < comp.nDValues(); i++)
    {
      addStdUnc(comp.StdUnc(i)*fact);
    }
  }else
  {
    uncertainty->setTypeToUnknown();
  }
// parameter with unc dealing (so value and Dvalue correspond)

  for(int i = 0; i < comp.nValues(); i++)
  {
    addValue(comp.value(i)*fact);
  }
}

void ParameterPhy::replace(const ParameterPhy &rhs)
{
  *this = rhs; 
  setNamePar(rhs.namePar());
}

//operators
ParameterPhy& ParameterPhy::operator= (const ParameterPhy &rhs)
{
/*!This is the assignement operator, not a replacement.
 * What it means is that the name, if existent, is not
 * replaced.
 *
 * Here we replace the values, uncertainty values and type and unit.
 */
  if(this == &rhs)return *this;
// first, the AtomicParameter
  parameter = rhs.getParameter();

// then, the UnitObject
  if(rhs.getUnitObject() != NULL)
  {
    copyUnitObject(rhs.getUnitObject());
  }else
  {
    if(unit == NULL)
    {
      unit = new Units;
    }else
    {
      unit->clear();
    }
  }

// finally, the UncObject
  if(uncertainty == NULL)uncertainty = new CoreUnc;
  if(rhs.getUncertaintyObject() != NULL)
  {
    *uncertainty = *(rhs.getUncertaintyObject());
  }else
  {
    uncertainty->setTypeToUnknown();
  }

  return *this;
}

ParameterPhy& ParameterPhy::operator+= (const ParameterPhy &rhs)
{
  const std::string method("ParameterPhy& ParameterPhy::operator+= (const ParameterPhy &)");
//// 1 - unit check
//// 1.a - test pointer
  if(rhs.getUnitObject() == NULL)
  {
    std::string errorStr = "It may be wrong to try to add physical parameters with units pointer NULL.\n";
    errorStr += "The adding part will be done, considering thus that this is a scalar, but there is a potential problem.";
    antiochWarning(method,errorStr);
  }
//// 1.b - unit homogeneity
  changeToSomeUnit(rhs.unitPar());

//// 2 - uncertainty
  modelUncertainty(rhs,&ParameterPhy::additiveModelUncertainty);

//// 3 - values
  parameter += rhs.getParameter();

  return *this;
}

ParameterPhy& ParameterPhy::operator-= (const ParameterPhy &rhs)
{
  *this += -rhs;
  return *this;
}

ParameterPhy& ParameterPhy::operator*= (const ParameterPhy &rhs)
{
  const std::string method("ParameterPhy& ParameterPhy::operator*= (const ParameterPhy &)");
//unit
  unit->add(rhs.unitPar());
//uncertainty
  modelUncertainty(rhs,&ParameterPhy::multiplyModelUncertainty);
//value
  parameter *= rhs.getParameter();

  return *this;
}

ParameterPhy & ParameterPhy::operator/= (const ParameterPhy &rhs)
{
//unit
  unit->substract(rhs.unitPar());
//unc
  modelUncertainty(rhs,&ParameterPhy::divideModelUncertainty);
//values
  parameter /= rhs.getParameter();

  return *this;
}

ParameterPhy & ParameterPhy::operator+= (double r)
{
  parameter += r;  
  return *this;
}

ParameterPhy & ParameterPhy::operator-= (double r)
{
  r *= -1.;
  *this += r;
  return *this;
}

ParameterPhy & ParameterPhy::operator*= (double r)
{
  parameter *= r;
  if(typeDPar() == CORE_UNCERTAINTY_TYPE_ABSOLUTE)
  {
     *uncertainty *= r;
  }
  return *this;
}

ParameterPhy & ParameterPhy::operator/= (double r)
{
  parameter /= r;
  if(typeDPar() == CORE_UNCERTAINTY_TYPE_ABSOLUTE)
  {
     *uncertainty /= r;
  }
  return *this;
}

ParameterPhy  ParameterPhy::operator+  (double r) const
{
  return (ParameterPhy(*this) += r);
}

ParameterPhy  ParameterPhy::operator-  (double r) const
{
  return(ParameterPhy(*this) -= r);
}

ParameterPhy  ParameterPhy::operator*  (double r) const
{
  return(ParameterPhy(*this) *= r);
}

ParameterPhy ParameterPhy::operator/  (double r) const
{
  return(ParameterPhy(*this) /= r);
}


ParameterPhy ParameterPhy::operator* (const ParameterPhy &rhs) const
{
  return (ParameterPhy(*this) *= rhs);
}

ParameterPhy ParameterPhy::operator/  (const ParameterPhy &rhs) const
{
  return (ParameterPhy(*this) /= rhs);
}

ParameterPhy ParameterPhy::operator+  (const ParameterPhy &rhs) const
{
  return (ParameterPhy(*this) += rhs);
}

ParameterPhy ParameterPhy::operator-  (const ParameterPhy &rhs) const
{
  return (ParameterPhy(*this) -= rhs);
}

ParameterPhy operator+ (const double &lhs, const ParameterPhy &rhs)
{
  ParameterPhy out(lhs);
  out += rhs;
  return out;
}
ParameterPhy operator- (const double &lhs, const ParameterPhy &rhs)
{
  ParameterPhy out(lhs);
  out -= rhs;
  return out;
}
ParameterPhy operator* (const double &lhs, const ParameterPhy &rhs)
{
  ParameterPhy out(lhs);
  out *= rhs;
  return out;
}
ParameterPhy operator/ (const double &lhs, const ParameterPhy &rhs)
{
  ParameterPhy out(lhs);
  out /= rhs;
  return out;
}
}
