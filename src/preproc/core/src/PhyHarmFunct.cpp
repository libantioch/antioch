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
#include "antioch/PhyHarmFunct.hpp"

namespace Antioch
{

PhyHarmFunct::PhyHarmFunct(const std::vector<PhyFunct> &functions, double xstep)
{
  abs = HarmonizeAbscissa(&functions,xstep);
  addOrd(functions);
}

PhyHarmFunct::PhyHarmFunct(const std::vector<PhyFunct> &functions, const std::string defcal, int xref, const std::string deftype)
{
  abs = HarmonizeAbscissa(&functions,xref,deftype);
  addOrd(functions,defcal);
}

PhyHarmFunct::PhyHarmFunct(const ParameterPhy &xgrid, const std::vector<PhyFunct> &functions, const std::string defcal):
  abs(xgrid)
{
  addOrd(functions,defcal);
}

void PhyHarmFunct::addOrd(const ParameterPhy &f)
{
   const std::string method("void PhyHarmFunct::addOrd(const ParameterPhy &)");
   if(abs.empty())antiochError(method,"Please set an abscissa before adding ordinate.");
   if(f.nValues() != abs.nValues())antiochError(method,"This ordinate (" 
                                                     + f.namePar() + 
                                                     ") does not have the same number of values than the abscissa ("
                                                     + abs.namePar() + ").");
  ord.push_back(f);
  hasOnlyZero.push_back(isZeroPP(f));
}

bool PhyHarmFunct::isZeroPP(const ParameterPhy &f)
{
  for(int i = 0; i < f.nValues(); i++)
  {
     if(f.value(i) != 0.)return false;
  }
  return true;
}

void PhyHarmFunct::addZeroOrd(const std::string &name,const std::string &units)
{
  ord.push_back(ParameterPhy(name,units,CORE_UNCERTAINTY_TYPE_NONE));
  for(int i = 0; i < abs.nValues(); i++)
  {
     ord.back().addValue(0.);
  }
  hasOnlyZero.push_back(true);
}

void PhyHarmFunct::setOrd(const PhyFunct &f,int iord, const std::string defcal)
{
  const std::string method("void PhyHarmFunct::setOrd(const PhyFunct &,int )");
  if(iord >= (int)ord.size())antiochError(method,"Requested ordinate over range of existing ordinates.");

  ord[iord] = setOrdToGrid(f,defcal);
  hasOnlyZero[iord] = isZeroPP(getOrd(iord));
}

bool PhyHarmFunct::abscissaHomogeneous(const std::vector<PhyFunct> *xs)
{
  std::string xUnit = xs->at(0).x.unitPar();
  for(unsigned int ix = 1; ix < xs->size(); ix++)
  {
     if(!xs->at(ix).x.isHomogeneous(xUnit))return false;
  }
  return true;
}

ParameterPhy PhyHarmFunct::HarmonizeAbscissa(const std::vector<PhyFunct> *xs, int xref, const std::string type)
{
  const std::string method("ParameterPhy PhyHarmFunct::HarmonizeAbscissa(const std::vector<PhyFunct> *, int, const std::string)");

//first, unit check
  std::string xUnit = xs->at(0).x.unitPar();
  if(!abscissaHomogeneous(xs))
      antiochError(method,"One abscissa (at least!) is not homogeneous to the others.\nReference used is in " + xUnit);

//second, chosen grid
  if(type == "abscissa") //if grid is given
  {
    if(xref >= (int)xs->size())antiochError(method,"Out-of-range reference abscissa requested");
    return xs->at(xref).x;
  }else //to be computed, (xmax - xmin)/nbins
  {
    ParameterPhy xgrid("grid",xUnit,CORE_UNCERTAINTY_TYPE_NONE);
    double xmin(1e303),xmax(-1e303);
    if(xref < 0)xref = 100; //default
    for(unsigned int ix = 0; ix < xs->size(); ix++)
    {
        for(int iv = 0; iv < xs->at(ix).x.nValues(); iv++)
        {
           if(xmax < xs->at(ix).x.value(iv))xmax = xs->at(ix).x.value(iv);
           if(xmin > xs->at(ix).x.value(iv))xmin = xs->at(ix).x.value(iv);
        }
    }
    for(int ig = 0; ig < xref; ig++)
    {
        xgrid.addValue(xmin + (double)ig * (xmax - xmin)/(double)(xref - 1));
    }
    return xgrid;
  }

}

ParameterPhy PhyHarmFunct::HarmonizeAbscissa(const std::vector<PhyFunct> *xs, double xstep)
{
  const std::string method("ParameterPhy PhyHarmFunct::HarmonizeAbscissa(const std::vector<PhyFunct> *, double )");

  std::string xUnit = xs->at(0).x.unitPar();
  if(!abscissaHomogeneous(xs))
      antiochError(method,"One abscissa (at least!) is not homogeneous to the others.\nReference used is in " + xUnit);

  ParameterPhy xgrid("grid",xUnit,CORE_UNCERTAINTY_TYPE_NONE);
  double xmin(1e303),xmax(-1e303);
  for(unsigned int ix = 0; ix < xs->size(); ix++)
  {
    for(int iv = 0; iv < xs->at(ix).x.nValues(); iv++)
    {
      if(xmax < xs->at(ix).x.value(iv))xmax = xs->at(ix).x.value(iv);
      if(xmin > xs->at(ix).x.value(iv))xmin = xs->at(ix).x.value(iv);
    }
  }

  for(double x = xmin; x <= xmax; x += xstep)
  {
    xgrid.addValue(x);
  }

  return xgrid; 
}

ParameterPhy PhyHarmFunct::setOrdToGrid(const PhyFunct &f, const std::string defcal)
{
  const std::string method("ParameterPhy PhyHarmFunct::setOrdToGrid(const PhyFunct &)");

  if(abs.empty())antiochError(method,"You MUST first fill the abscissa before finding the ordinate at those missing abscissa values.");

  if(defcal == "interpolation")
  {
//value by linear integration, dvalue by conservative hypothesis: largest
    return interpolationOnGrid(f);
  }else
  {
    return integrationOnGrid(f);
  }

}

ParameterPhy PhyHarmFunct::integrationOnGrid(const PhyFunct &f)
{
  ParameterPhy out(f.y.namePar(),f.y.unitPar(),f.y.typeDPar());
  bool computeDVal(!f.y.noUnc() && !f.y.unknownUnc());
  double coef(f.x.getFactorToSomeUnit(abs.unitPar()));
  int iy(0);
  while(coef*f.x.value(iy) < abs.value(0))iy++;
  if(coef*f.x.value(iy) > abs.value(1) && iy > 0)iy--;
  for(int iv = 0; iv < abs.nValues() - 1; iv++)
  {
     double mean(0.),dmean(0.);
     int np(0);
     while((coef*f.x.value(iy) < abs.value(iv + 1)) && (iy < f.x.nValues() - 1))
     {
        mean  += f.y.value(iy);
        dmean += f.y.variance(iy);
        iy++;
        np++;
     }
     if(np > 0)//mean if data points in interval
     {
       mean  /= (double)np;
       dmean /= (double)np;
     }else//interpolation if not, the iy and next to iy points.
     {
       double a,b;
       int iiy(iy + 1);
       if(iiy >= f.x.nValues() || abs.value(iv) < coef*f.x.value(iy))iiy = iy - 1;
       interpolation(a,b,coef*f.x.value(iy),f.y.value(iy),coef*f.x.value(iiy),f.y.value(iiy));
       mean = a * abs.value(iv) + b;
       dmean = f.y.variance(iy);
       if(f.y.variance(iy) < f.y.variance(iiy))dmean = f.y.variance(iiy);
     }
     out.addValue(mean);
     if(computeDVal)out.addVariance(dmean);
  }
//last value
  double a,b;
  interpolation(a,b,coef*f.x.value(f.x.nValues()- 2),f.y.value(f.y.nValues()- 2),coef*f.x.value(f.x.nValues()- 1),f.y.value(f.y.nValues()- 1));
  out.addValue(a*abs.value(abs.nValues() - 1) + b);

  return out;
}

ParameterPhy PhyHarmFunct::interpolationOnGrid(const PhyFunct &f)
{
  const std::string method("ParameterPhy PhyHarmFunct::interpolationOnGrid(const PhyFunct &)");

  ParameterPhy out(f.y.namePar(),f.y.unitPar(),f.y.typeDPar());
  bool computeDVal(!f.y.noUnc() && !f.y.unknownUnc());

  double coef(f.x.getFactorToSomeUnit(abs.unitPar()));

  double xmin(1e303),xmax(-1e303);
  for(int iv = 0; iv < f.x.nValues(); iv++)
  {
    if(xmin > f.x.value(iv))xmin = f.x.value(iv);
    if(xmax < f.x.value(iv))xmax = f.x.value(iv);
  }

  for(int iv = 0; iv < abs.nValues(); iv++)
  {
    if(abs.value(iv) < xmin*coef || abs.value(iv) > xmax*coef)
    {
      out.addValue(0.);
      if(computeDVal)out.addDValue(0.);
      continue;
    }
    double xinf(1e303),xsup(1e303);
    int fvinf(-1),fvsup(-1);
    for(int fv = 0; fv < f.x.nValues(); fv++)
    {
      if(f.x.value(fv)*coef == abs.value(iv))
      {
        fvinf = fv;
        fvsup = fv;
        break;
      }else if(f.x.value(fv)*coef < abs.value(iv))
      {
        if(xinf > (abs.value(iv) - f.x.value(fv)*coef))
        {
           xinf = abs.value(iv) - f.x.value(fv)*coef;
           fvinf = fv;
        }
      }else
      {
        if(xsup > (f.x.value(fv)*coef - abs.value(iv)))
        {
           xsup = f.x.value(fv)*coef - abs.value(iv);
           fvsup = fv;
        }
      }
    }
    if(fvinf == -1 || fvsup == -1)
    {
      antiochError(method,"Did not find an inferior or superior value");
    }else if(fvinf == fvsup)
    {
      out.addValue(f.y.value(fvinf));
      if(computeDVal)out.addDValue(f.y.Dvalue(fvinf));
    }else
    {
//linear extrapolation
      double a,b;
      interpolation(a,b,f.x.value(fvinf),f.y.value(fvinf),f.x.value(fvsup),f.y.value(fvsup));
      out.addValue(a*abs.value(iv)/coef + b);

//the biggest
      if(computeDVal)
      {
        double dval(f.y.Dvalue(fvinf));
        if(dval < f.y.Dvalue(fvsup))dval = f.y.Dvalue(fvsup);
        out.addDValue(dval);
      }
    }
  }

  return out;
}

void PhyHarmFunct::interpolation(double &a, double &b, double x1, double y1, double x2, double y2)
{
  a = (y2 - y1)/(x2 - x1);
  b = y2 - (y2 - y1)/(x2 - x1) * x2;
}

void PhyHarmFunct::addOrd(const PhyFunct &f, const std::string defcal)
{
  ord.push_back(setOrdToGrid(f,defcal));
  hasOnlyZero.push_back(isZeroPP(ord.back()));
}

void PhyHarmFunct::addOrd(const std::vector<PhyFunct> &f, const std::string defcal)
{
  for(unsigned int i = 0; i < f.size(); i++)
  {
     addOrd(f[i],defcal);
  }
}

void PhyHarmFunct::addOrd(const PhyHarmFunct &f, int iord, const std::string defcal)
{
  PhyFunct tmp;
  tmp.x = f.getAbs();
  tmp.y = f.getOrd(iord);

  addOrd(tmp,defcal);
}

int PhyHarmFunct::indexByName(const std::string &nameOrd)
{
  const std::string method("int PhyHarmFunct::indexByName(const std::string &)");
  for(unsigned int i = 0; i < ord.size(); i++)
  {
     if(ord[i].namePar() == nameOrd)return i;
  }

  antiochError(method,"Did not find ordinate of name "+ nameOrd);
  return -1;
}

ParameterPhy PhyHarmFunct::getOrd(const std::string &nameOrd) const
{
  const std::string method("ParameterPhy PhyHarmFunct::getOrd(const std::string &)");
  for(unsigned int i = 0; i < ord.size(); i++)
  {
     if(ord[i].namePar() == nameOrd)return ord[i];
  }

  antiochError(method,"Did not find ordinate of name "+ nameOrd);
  return ParameterPhy();
}

void PhyHarmFunct::showAll(std::ostream &out) const
{
  out << "#PhyHarmFunct" << std::endl;
  out << "#name: " << name << std::endl;
  out << "Abscissa" << std::endl;
  abs.showAll(out);
  out << "Ordinates" << std::endl;
  for(int i = 0; i < nOrdinates(); i++)
  {
    ord[i].showAll(out);
  }
}

PhyHarmFunct & PhyHarmFunct::operator= (const PhyHarmFunct &rhs)
{
  setName(rhs.getName());
  setAbs(rhs.getAbs());
  ord = rhs.getAllOrd();

  return *this;
}


bool PhyHarmFunct::operator== (const PhyHarmFunct &rhs)
{
  if(getName() != rhs.getName())return false;
  if(getAbs() != rhs.getAbs())return false;
  if(nOrdinates() != rhs.nOrdinates())return false;
  for(int i = 0; i < nOrdinates(); i++)
  {
    if(getOrd(i) != rhs.getOrd(i))return false;
  }

  return true;
}
}
