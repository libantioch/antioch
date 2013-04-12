//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

//Antioch
#include "antioch/phymath.hpp"

//C++
#include <cmath>

namespace Antioch{
/********************************************
 *
 * First the generic methods
 *
 ********************************************/

double checkFunctionUnit(const std::string unitRef, 
                         const std::string methodCalling, 
                         const std::string errmessage,
                         const ParameterPhy &rhs)
{
  double coef(1.);
  if(rhs.unitPar() != unitRef)
  {
     Units *tmp = new Units(rhs.unitPar());
     if(!tmp->isHomogeneous(unitRef))
     {
        delete tmp;
        antiochError(methodCalling,errmessage);
        return 0.;
     }
     coef = tmp->FactorToSomeUnit(unitRef);
     delete tmp;
  }

  return coef;
}

void calculateFunctionValues(ParameterPhy &y, 
                             const ParameterPhy &x, 
                             double coef, 
                             double function(double ))
{
  for(int i = 0; i < x.nValues(); i++)
  {
    y.addValue(function(x.value(i)*coef));
  }
  
}

void calculateFunctionDValues(ParameterPhy &y, 
                             const ParameterPhy &x, 
                             double coef, 
                             double function(double ))
{
  if(!x.isUncDef())
  {
     y.setTypeDPar(x.typeDPar());
     return;
  }

  y.setTypeDPar(CORE_UNCERTAINTY_TYPE_ABSOLUTE);
  if(x.typeDPar() == CORE_UNCERTAINTY_TYPE_RELATIVE)
  {
    for(int i = 0; i < x.nDValues(); i++)
    {
       y.addDValue(function(x.value(i)*coef) * function(x.value(i)*coef) * x.Dvalue(i)*x.value(i)*coef * x.Dvalue(i)*x.value(i)*coef);
    }
  }else{
    for(int i = 0; i < x.nDValues(); i++)
    {
       y.addDValue(function(x.value(i)*coef) * function(x.value(i)*coef) * x.variance(i)*coef*coef);
    }
  }
}

/***********************************************
 *
 * Now the double f(double x) needed methods
 *
 ***********************************************/
double dtandx(double x)
{
  return 1./(std::cos(x)*std::cos(x));
}

double dacosdx(double x)
{
  return -1./(std::sqrt(1. - x*x));
}

double dasindx(double x)
{
  return 1./(std::sqrt(1. - x*x));
}

double datandx(double x)
{
  return 1./(1. + x*x);
}

double dlogdx(double x)
{
  return 1./x;
}

double dlog10dx(double x)
{
  return 1./(std::log(10.)*x);
}

double erf(double z)
{
  if(z < 0.)z = -z;
  double p(0.3275911);
  double t(1./(1. + p*z));
  double a1(0.254829592);
  double a2(-0.284496736);
  double a3(1.421413741);
  double a4(-1.453152027);
  double a5(1.061405429);

  return (1. - (a1*t + a2*t*t + a3*t*t*t + a4*t*t*t*t + a5*t*t*t*t*t) * std::exp(-z*z));
}

double derfdx(double x)
{
  return (2./std::sqrt(3.141592653589793) * std::exp(-x*x));
}

double dsqrtdx(double x)
{
  return 1./(2.*std::sqrt(x));
}

/***********************************************
 *
 * Now the ParameterPhy f(ParameterPhy x) needed methods
 *
 ***********************************************/

ParameterPhy cos(const ParameterPhy &rhs)
{
  const std::string method("ParameterPhy cos(const ParameterPhy &)");
  const std::string errmessage("Parameter is not an angle");
  const std::string unitRef("rad");

//unit management
  double coef = checkFunctionUnit(unitRef,method,errmessage,rhs);

  ParameterPhy out;
  out.setUnitPar("");

//value calculations
  calculateFunctionValues(out,rhs,coef,std::cos);

//unc management and calculations if necessary
  calculateFunctionDValues(out,rhs,coef,std::sin);


  return out;
}

ParameterPhy sin(const ParameterPhy &rhs)
{
  const std::string method("ParameterPhy sin(const ParameterPhy &)");
  const std::string errmessage("Parameter is not an angle");
  const std::string unitRef("rad");

//unit management
  double coef = checkFunctionUnit(unitRef,method,errmessage,rhs);

  ParameterPhy out;
  out.setUnitPar("");

//value calculations
  calculateFunctionValues(out,rhs,coef,std::sin);

//unc management and calculations if necessary
  calculateFunctionDValues(out,rhs,coef,std::cos);


  return out;
}

ParameterPhy tan(const ParameterPhy &rhs)
{
  const std::string method("ParameterPhy tan(const ParameterPhy &)");
  const std::string errmessage("Parameter is not an angle");
  const std::string unitRef("rad");

//unit management
  double coef = checkFunctionUnit(unitRef,method,errmessage,rhs);

  ParameterPhy out;
  out.setUnitPar("");

//calculations
  calculateFunctionValues(out,rhs,coef,std::tan);

//unc management and calculations if necessary
  calculateFunctionDValues(out,rhs,coef,dtandx);


  return out;
}

ParameterPhy acos(const ParameterPhy &rhs)
{
  const std::string method("ParameterPhy acos(const ParameterPhy &)");
  const std::string errmessage("Parameter has unit " + rhs.unitPar());
  const std::string unitRef("");

//unit management
  double coef = checkFunctionUnit(unitRef,method,errmessage,rhs);

  ParameterPhy out;
  out.setUnitPar("rad");

//calculations
  calculateFunctionValues(out,rhs,coef,std::acos);


//unc management and calculations if necessary
  calculateFunctionDValues(out,rhs,coef,dacosdx);


  return out;
}

ParameterPhy asin(const ParameterPhy &rhs)
{
  const std::string method("ParameterPhy asin(const ParameterPhy &)");
  const std::string errmessage("Parameter has unit " + rhs.unitPar());
  const std::string unitRef("");

//unit management
  double coef = checkFunctionUnit(unitRef,method,errmessage,rhs);

  ParameterPhy out;
  out.setUnitPar("rad");

//calculations
  calculateFunctionValues(out,rhs,coef,std::asin);

//unc management and calculations if necessary
  calculateFunctionDValues(out,rhs,coef,dasindx);


  return out;
}

ParameterPhy atan(const ParameterPhy &rhs)
{
  const std::string method("ParameterPhy atan(const ParameterPhy &)");
  const std::string errmessage("Parameter has unit " + rhs.unitPar());
  const std::string unitRef("");

//unit management
  double coef = checkFunctionUnit(unitRef,method,errmessage,rhs);

  ParameterPhy out;
  out.setUnitPar("rad");

//calculations
  calculateFunctionValues(out,rhs,coef,std::atan);

//unc management and calculations if necessary
  calculateFunctionDValues(out,rhs,coef,datandx);


  return out;
}

ParameterPhy exp(const ParameterPhy &rhs)
{
  const std::string method("ParameterPhy exp(const ParameterPhy &)");
  const std::string errmessage("Parameter has unit " + rhs.unitPar());
  const std::string unitRef("");

//unit management
  double coef = checkFunctionUnit(unitRef,method,errmessage,rhs);

  ParameterPhy out;
  out.setUnitPar("");

//calculations
  calculateFunctionValues(out,rhs,coef,std::exp);

//unc management and calculations if necessary
  calculateFunctionDValues(out,rhs,coef,std::exp);


  return out;
}

ParameterPhy log(const ParameterPhy &rhs)
{
  const std::string method("ParameterPhy log(const ParameterPhy &)");
  const std::string errmessage("Parameter has unit " + rhs.unitPar());
  const std::string unitRef("");

//unit management
  double coef = checkFunctionUnit(unitRef,method,errmessage,rhs);

  ParameterPhy out;
  out.setUnitPar("");

//calculations
  calculateFunctionValues(out,rhs,coef,std::log);

//unc management and calculations if necessary
  calculateFunctionDValues(out,rhs,coef,dlogdx);


  return out;
}

ParameterPhy log10(const ParameterPhy& rhs)
{
  const std::string method("ParameterPhy log10(const ParameterPhy &)");
  const std::string errmessage("Parameter has unit " + rhs.unitPar());
  const std::string unitRef("");

//unit management
  double coef = checkFunctionUnit(unitRef,method,errmessage,rhs);

  ParameterPhy out;
  out.setUnitPar("");

//calculations
  calculateFunctionValues(out,rhs,coef,std::log10);

//unc management and calculations if necessary
  calculateFunctionDValues(out,rhs,coef,dlog10dx);


  return out;
}


ParameterPhy fabs(const ParameterPhy &rhs)
{
  ParameterPhy out;
//unit management, same unit
  out.setUnitPar(rhs.unitPar());
  double coef(1.);

//calculations
  calculateFunctionValues(out,rhs,coef,std::fabs);

//unc management and calculations if necessary, no need to go into trouble, just copying
  out.setTypeDPar(rhs.typeDPar());
  for(int i = 0; i < rhs.nDValues(); i++)
  {
     out.addDValue(rhs.Dvalue(i));
  }


  return out;
}

ParameterPhy erf(const ParameterPhy& rhs)
{
  const std::string method("ParameterPhy erf(const ParameterPhy &)");
  const std::string errmessage("Parameter has unit " + rhs.unitPar());
  const std::string unitRef("");

//unit management
  double coef = checkFunctionUnit(unitRef,method,errmessage,rhs);

  ParameterPhy out;
  out.setUnitPar("");

//calculations
  calculateFunctionValues(out,rhs,coef,erf);

//unc management and calculations if necessary
  calculateFunctionDValues(out,rhs,coef,derfdx);


  return out;
}

ParameterPhy sqrt(const ParameterPhy& rhs)
{
  const std::string method("ParameterPhy sqrt(const ParameterPhy &)");

//unit management
  Units *unit = new Units(rhs.unitPar());
  unit->root(2);

  ParameterPhy out;
  out.setUnitPar(unit->getSymbol());
  delete unit;

//calculations
  calculateFunctionValues(out,rhs,1.,std::sqrt);

//unc management and calculations if necessary
  calculateFunctionDValues(out,rhs,1.,dsqrtdx);


  return out;
}


ParameterPhy pow(const ParameterPhy& rhs, int p)
{
  ParameterPhy out;
//unit management
  Units *outUnit = new Units(rhs.unitPar());
  *outUnit *= p;
  out.setUnitObject(outUnit);

//values calculations
  for(int i = 0; i < rhs.nValues(); i++)
  {
    out.addValue(std::pow(rhs.value(i),p));
  }

//uncertainty management
  if(rhs.noUnc() || rhs.unknownUnc())
  {
      out.setTypeDPar(rhs.typeDPar());
      return out;
  }

//uncertainty calculations
  out.setTypeDPar(CORE_UNCERTAINTY_TYPE_ABSOLUTE);
  if(rhs.typeDPar() == CORE_UNCERTAINTY_TYPE_RELATIVE)
  {
    for(int i = 0; i < rhs.nValues(); i++)
    {
      out.addDValue(rhs.Dvalue(i)*rhs.value(i)*rhs.Dvalue(i)*rhs.value(i)*(double)p*std::pow(rhs.value(i),p-1)*(double)p*std::pow(rhs.value(i),p-1));
    }
  }else
  {
    for(int i = 0; i < rhs.nValues(); i++)
    {
      out.addDValue(rhs.Dvalue(i)*(double)p*std::pow(rhs.value(i),p-1)*(double)p*std::pow(rhs.value(i),p-1));
    }
  }

  return out;
}

ParameterPhy pow(const ParameterPhy &rhs, double p)
{
  const std::string method("ParameterPhy pow(const ParameterPhy &, double)");
  const std::string errmessage("Parameter has unit " + rhs.unitPar());
  const std::string unitRef("");

//unit management
  double coef = checkFunctionUnit(unitRef,method,errmessage,rhs);
  ParameterPhy out;
  out.setUnitPar("");

//calculations
  for(int i = 0; i < rhs.nValues(); i++)
  {
    out.addValue(std::pow(rhs.value(i)*coef,p));
  }

//uncertainty management
  if(rhs.noUnc() || rhs.unknownUnc())
  {
      out.setTypeDPar(rhs.typeDPar());
      return out;
  }

//uncertainty calculations
  out.setTypeDPar(CORE_UNCERTAINTY_TYPE_ABSOLUTE);
  if(rhs.typeDPar() == CORE_UNCERTAINTY_TYPE_RELATIVE)
  {
    for(int i = 0; i < rhs.nValues(); i++)
    {
      out.addDValue(rhs.Dvalue(i)*rhs.value(i)*coef*rhs.Dvalue(i)*rhs.value(i)*coef*p*std::pow(rhs.value(i)*coef,p-1.)*p*std::pow(rhs.value(i)*coef,p-1.));
    }
  }else
  {
    for(int i = 0; i < rhs.nValues(); i++)
    {
      out.addDValue(rhs.Dvalue(i)*coef*coef*p*std::pow(rhs.value(i)*coef,p-1.)*p*std::pow(rhs.value(i)*coef,p-1.));
    }
  }

  return out;
}

ParameterPhy pow(double b, const ParameterPhy &rhs)
{
  const std::string method("ParameterPhy pow(double, const ParameterPhy &)");
  const std::string errmessage("Parameter has unit " + rhs.unitPar());
  const std::string unitRef("");

//unit management
  double coef = checkFunctionUnit(unitRef,method,errmessage,rhs);
  ParameterPhy out;
  out.setUnitPar("");

//calculations
  for(int i = 0; i < rhs.nValues(); i++)
  {
    out.addValue(std::pow(b,rhs.value(i)*coef));
  }

//uncertainty management
  if(rhs.noUnc() || rhs.unknownUnc())
  {
      out.setTypeDPar(rhs.typeDPar());
      return out;
  }

//uncertainty calculations
  out.setTypeDPar(CORE_UNCERTAINTY_TYPE_ABSOLUTE);
  if(rhs.typeDPar() == CORE_UNCERTAINTY_TYPE_RELATIVE)
  {
    for(int i = 0; i < rhs.nValues(); i++)
    {
      out.addDValue(std::log(b)*std::pow(b,rhs.value(i)*coef)*std::log(b)*std::pow(b,rhs.value(i)*coef)*rhs.Dvalue(i)*coef*rhs.Dvalue(i)*coef*rhs.value(i)*coef);
    }
  }else
  {
    for(int i = 0; i < rhs.nValues(); i++)
    {
      out.addDValue(std::log(b)*std::pow(b,rhs.value(i)*coef)*std::log(b)*std::pow(b,rhs.value(i)*coef)*rhs.Dvalue(i)*coef*rhs.Dvalue(i)*coef);
    }
  }

  return out;
}

ParameterPhy pow(const ParameterPhy &base, const ParameterPhy &p)
{
  const std::string method("ParameterPhy pow(const ParameterPhy &, const ParameterPhy &)");
  const std::string baseErrmessage("The base of this power operation is not unitless: " + base.unitPar());
  const std::string powerErrmessage("The powe of this power operation is not unitless: " + p.unitPar());
  const std::string unitRef("");

//unit management
  double coefbase  = checkFunctionUnit(unitRef,method,baseErrmessage,base);
  double coefpower = checkFunctionUnit(unitRef,method,powerErrmessage,p);
  ParameterPhy out;
  out.setUnitPar("");
//unit management

//calculations
  if(base.nValues() != 1 && p.nValues() != 1 && base.nValues() != p.nValues())
  {
    antiochError(method,"Base and exposant parameters have different numbers of value, and no one has only one.");
    return out;
  }
  int nval;
  (p.nValues() > base.nValues())?nval = p.nValues():nval = base.nValues();
  for(int i = 0; i < nval; i++)
  {
     int baseval(i),pval(i);
     if(base.nValues() == 1)baseval=0;
     if(p.nValues() == 1)pval=0;
     out.addValue(std::pow(base.value(baseval)*coefbase,p.value(pval)*coefpower));
  }

//uncertainty management
  if(base.noUnc() || base.unknownUnc())
  {
     out.setTypeDPar(base.typeDPar());
     return out;
  }else if(p.noUnc() || p.unknownUnc())
  {
     out.setTypeDPar(p.typeDPar());
     return out;
  }

//uncertainty calculations
  out.setTypeDPar(CORE_UNCERTAINTY_TYPE_ABSOLUTE);
  if(base.typeDPar() == CORE_UNCERTAINTY_TYPE_RELATIVE && p.typeDPar() == CORE_UNCERTAINTY_TYPE_RELATIVE)
  {
    for(int i = 0; i < nval; i++)
    {
     int baseval(i),pval(i);
     if(base.nValues() == 1)baseval=0;
     if(p.nValues() == 1)pval=0;
// f = x^y
     out.addDValue(// dx^2 * y^2 * (x^(y-1))^2
                base.Dvalue(baseval)*base.value(baseval)*coefbase*base.Dvalue(baseval)*base.value(baseval)*coefbase
                   *p.value(pval)*coefpower*p.value(pval)*coefpower*std::pow(base.value(baseval)*coefbase,2.*(p.value(pval)*coefpower - 1.))
                + //dy^2 * (ln(x)*x^y)^2
                std::log(base.value(baseval)*coefbase)*std::log(base.value(baseval)*coefbase)*std::pow(base.value(baseval)*coefbase,p.value(pval)*coefpower)
                        *std::pow(base.value(baseval)*coefbase,p.value(pval)*coefpower)
                        *p.Dvalue(pval)*p.value(pval)*coefpower*p.Dvalue(pval)*p.value(pval)*coefpower
                  );
    }
  }else if(base.typeDPar() == CORE_UNCERTAINTY_TYPE_RELATIVE)
  {
    for(int i = 0; i < nval; i++)
    {
     int baseval(i),pval(i);
     if(base.nValues() == 1)baseval=0;
     if(p.nValues() == 1)pval=0;
// f = x^y
     out.addDValue(// dx^2 * y^2 * (x^(y-1))^2
                base.Dvalue(baseval)*base.value(baseval)*coefbase*base.Dvalue(baseval)*base.value(baseval)*coefbase
                   *p.value(pval)*coefpower*p.value(pval)*coefpower*std::pow(base.value(baseval)*coefbase,2.*(p.value(pval)*coefpower - 1.))
                + //dy^2 * (ln(x)*x^y)^2
                std::log(base.value(baseval)*coefbase)*std::log(base.value(baseval)*coefbase)*std::pow(base.value(baseval)*coefbase,p.value(pval)*coefpower)
                        *std::pow(base.value(baseval)*coefbase,p.value(pval)*coefpower)
                        *p.Dvalue(pval)*coefpower*coefpower
                  );
    }
  }else if(p.typeDPar() == CORE_UNCERTAINTY_TYPE_RELATIVE)
  {
    for(int i = 0; i < nval; i++)
    {
     int baseval(i),pval(i);
     if(base.nValues() == 1)baseval=0;
     if(p.nValues() == 1)pval=0;
// f = x^y
     out.addDValue(// dx^2 * y^2 * (x^(y-1))^2
                base.Dvalue(baseval)*coefbase*coefbase
                   *p.value(pval)*coefpower*p.value(pval)*coefpower*std::pow(base.value(baseval)*coefbase,2.*(p.value(pval)*coefpower - 1.))
                + //dy^2 * (ln(x)*x^y)^2
                std::log(base.value(baseval)*coefbase)*std::log(base.value(baseval)*coefbase)*std::pow(base.value(baseval)*coefbase,p.value(pval)*coefpower)
                        *std::pow(base.value(baseval)*coefbase,p.value(pval)*coefpower)
                        *p.Dvalue(pval)*p.value(pval)*coefpower*p.Dvalue(pval)*p.value(pval)*coefpower
                  );
    }
  }else
  {
    for(int i = 0; i < nval; i++)
    {
     int baseval(i),pval(i);
     if(base.nValues() == 1)baseval=0;
     if(p.nValues() == 1)pval=0;
// f = x^y
     out.addDValue(// dx^2 * y^2 * (x^(y-1))^2
                base.Dvalue(baseval)*coefbase*coefbase
                   *p.value(pval)*coefpower*p.value(pval)*coefpower*std::pow(base.value(baseval)*coefbase,2.*(p.value(pval)*coefpower - 1.))
                + //dy^2 * (ln(x)*x^y)^2
                std::log(base.value(baseval)*coefbase)*std::log(base.value(baseval)*coefbase)*std::pow(base.value(baseval)*coefbase,p.value(pval)*coefpower)
                        *std::pow(base.value(baseval)*coefbase,p.value(pval)*coefpower)
                        *p.Dvalue(pval)*coefpower*coefpower
                  );
    }
  }


  return out;
}
}
