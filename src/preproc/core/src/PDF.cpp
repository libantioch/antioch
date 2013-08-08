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

//Antioch
#include "antioch/PDF.hpp"
#include "antioch/Error.hpp"

//C++
#include <sstream>
#include <cmath>

namespace Antioch{

Pdf *objectPDF(const std::string &type,const std::vector<double> &parameters, const std::vector<std::vector<double> > &mat)
{
  const std::string method("Pdf *objectPDF(const std::string &,std::vector<double> & (= std::vector<double>()) )");

  Pdf *ptr = NULL;

  if(type == KNOWN_PDF[0])
  {
     ptr = new Norm(parameters,mat);
  }else if(type == KNOWN_PDF[1])
  {
     ptr = new NorT(parameters,mat);
  }else if(type == KNOWN_PDF[2])
  {
     ptr = new Delta(parameters);
  }else if(type == KNOWN_PDF[3])
  {
     ptr = new Unif(parameters,mat);
  }else if(type == KNOWN_PDF[4])
  {
     ptr = new LogN(parameters,mat);
  }else if(type == KNOWN_PDF[5])
  {
     ptr = new LogU(parameters,mat);
  }else if(type == KNOWN_PDF[6])
  {
     ptr = new Diri(parameters);
  }else if(type == KNOWN_PDF[7])
  {
     ptr = new DirG(parameters);
  }else if(type == KNOWN_PDF[8])
  {
     ptr = new DiUT(parameters);
  }else if(type == KNOWN_PDF[9])
  {
     int input = (int)parameters.front();
     ptr = new DiUn(input);
  }else if(type == KNOWN_PDF[10])
  {
     std::vector<int> parameterInt;
     for(unsigned int i = 0; i < parameters.size(); i++)
     {
        parameterInt.push_back((int)parameters[i]);
     }
     ptr = new DiOr(parameterInt);
  }else
  {
     antiochWarning(method,"Undefined pdf for this parameter. " + type + " is not (yet) supported. NULL pointer returned.");
  }

  return ptr;
}


double Marginal::getParameter(int n) const
{
  const std::string method("double Marginal::getParameter(int ) const");
  if(n >= (int)parameters.size())
        antiochError(method,"Out-of-range parameter asked.");

  return parameters[n];
}


Marginal Pdf::marginalFromIndex(int ind)
{
  const std::string method("Marginal Pdf::marginalFromIndex(int )");

  if(ind >= (int)marginal.size())
  {
     antiochError(method,"Requested member index is beyond the numbers of stored marginal(s)");
     return marginal[0];//keeps compiler happy
  }else
  {
     return marginal[ind];
  }

}

Marginal Pdf::marginalFromName(const std::string &name) const
{
  const std::string method("Marginal Pdf::marginalFromName(const std::string &)");

  for(unsigned int i = 0; i < marginal.size(); i++)
  {
     if(name == getNameMar(i))return marginal[i];
  }

  antiochError(method,"Have not found marginal " + name + " in the current pdf.");

//So the compiler doesn't complain
  return marginal[0];
}

const std::vector<double> Pdf::getMeanVec() const
{

  std::vector<double> means;
  for(unsigned int i = 0; i < marginal.size(); i++)
  {
     means.push_back(marginal[i].getMean());
  }

  return means;
}

bool Pdf::isItReady(std::string &error)
{
  if(namePdf.find("Non initialized") != std::string::npos)
  {
    error = "It has not been initialized (namePdf says so).";
    return false;
  }
  if(marginal.size() == 0)
  {
    error = "No marginal.";
    return false;
  }
  for(unsigned int i = 0; i < marginal.size(); i++)
  {
     if(marginal[i].getMean() == -1e303)
     {
       error = "A marginal has no mean.";
       return false;
     }
     if(marginal[i].getMode() == -1e303)
     {
       error = "A marginal has no mode.";
       return false;
     }
     if(marginal[i].getMedian() == -1e303)
     {
       error = "A marginal has no median.";
       return false;
     }
     if(marginal[i].getStdDev() == -1e303)
     {
       error = "A marginal has no standard deviation.";
       return false;
     }
     if(marginal[i].getUpLimit() == -1e303)
     {
       error = "A marginal has no upper bound.";
       return false;
     }
     if(marginal[i].getDownLimit() == -1e303)
     {
       error = "A marginal has no lower bound.";
       return false;
     }
  }

  return true;
}

void Pdf::isItGood()
{
  std::string error;
  const std::string method("void Pdf::isItGood()");
  if(!isItReady(error))
  {
    std::string errStr = "This pdf (" + namePdf + ") was clearly not initialized properly.";
    if(marginal.size() == 0)errStr += " " + error;
    antiochError(method,errStr);
  }
}

void Pdf::setPDF(const std::vector<double>&data, const std::vector<std::vector<double> > &mat)
{
  const std::string method("void Pdf::setPDF(const std::vector<double>&, const std::vector<std::vector<double> > & (= std::vector<std::vector<double> >()) );");
  if(namePdf.find("Non initialized ") != std::string::npos)namePdf.erase(0,16);

  marginal.clear();
  addMarginalFromData(data);
  
  if(!mat.empty())//if given
  {
    Cov = mat;
  }else//else it's uncorrelated
  {
    Cov.resize(marginal.size());
    for(unsigned int i = 0; i < marginal.size(); i++)
    {
       Cov[i].resize(marginal.size());
       Cov[i][i] = 1.;
    }
  }
}

void Pdf::addMarginalFromData(const std::vector<double>&data)
{
  const std::string method("void Pdf::addMarginalFromData(const std::vector<double>&)");
  antiochError(method,"This pdf (" + namePdf + ") is not a pdf to be added marginal(s) to, you should not be here.");
}

void Pdf::testNumberOfParameters(unsigned int nP)
{
  const std::string method("void Pdf::testNumberOfParameters(unsigned int)");
  antiochError(method,"This pdf ("+ namePdf +") is not a pdf to be tested, you should not be here.");
}

void Pdf::initializePdf()
{
  const std::string method("void Pdf::initializePdf()");
  antiochError(method,"This pdf ("+ namePdf +") is not a pdf to be initialized, you should not be here.");
}

Norm::Norm(double m, double stdUnc, const std::string &nP,const std::string &id)
        :Pdf(nP,id)
{
  std::vector<double> data;
  data.push_back(m);
  data.push_back(stdUnc);
  setPDF(data);
}

void Norm::testNumberOfParameters(unsigned int nP)
{
  const std::string method("void Norm::testNumberOfParameters(unsigned int )");
  if(nP % 2 != 0)
        antiochError(method,"I need an even number of parameters for a normal pdf (mean and std dev for all marginals).");
}

void Norm::addMarginalFromData(const std::vector<double>&data)
{
  testNumberOfParameters(data.size());
  
  for(unsigned int i = 0; i < data.size(); i += 2)
  {
    marginal.push_back(Marginal(data[i], //mean
                                data[i], //mode
                                data[i], //median
                                data[i+1], //std dev
                                data[i] + 4.*data[i+1], //up
                                data[i] - 4.*data[i+1], //low
                                IdPdf) //same id
                      );
  }
}

NorT::NorT(double m, double stdUnc, double up, double down):
  Norm(m,stdUnc,"Normal truncated pdf","nort")
{
  setUpLimitMar(up);
  setDownLimitMar(down);
}

NorT::NorT(const std::vector<double> &norm, const std::vector<double> &limits, const std::vector<std::vector<double> >&mat, const std::string &nP, const std::string &id):
  Norm(norm,mat,nP,id)
{
  for(unsigned int i = 0; i < limits.size() && i < 2*(unsigned int)getNMarginal(); i += 2)
  {
    setUpLimitMar(limits[i],i);
    setDownLimitMar(limits[i+1],i);
  }
}



void NorT::testNumberOfParameters(unsigned int nP)
{
  const std::string method("void NorT::testNumberOfParameters(unsigned int )");
  if(nP % 4 != 0)
  {
    std::string errStr = "I need an number of parameters that is a multiple of 4";
    errStr += " for a normal truncated pdf (mean, std dev, upper bound and lower bound for all marginals).";
    antiochError(method,errStr);
  }
}

void NorT::addMarginalFromData(const std::vector<double>&data)
{
  
  testNumberOfParameters(data.size());
  for(unsigned int i = 0; i < data.size(); i += 4)
  {
    marginal.push_back(Marginal(data[i], //mean
                                data[i], //mode
                                data[i], //median
                                data[i+1], //std dev
                                data[i + 2], //up
                                data[i + 3], //low
                                IdPdf) //same id
                      );
  }
}

Unif::Unif(double up,double down,const std::string &nP, const std::string &id):
  Pdf(nP,id)
{
  std::vector<double> data;
  data.push_back(up);
  data.push_back(down);
  setPDF(data);
}

void Unif::testNumberOfParameters(unsigned int nP)
{
  const std::string method("void Unif::testNumberOfParameters(unsigned int )");
  if(nP % 2 != 0)
        antiochError(method,"I need an even number of parameters for a uniform pdf (upper and lower bounds for all marginals).");
}

void Unif::addMarginalFromData(const std::vector<double>&data)
{
  testNumberOfParameters(data.size());
  
  for(unsigned int i = 0; i < data.size(); i += 2)
  {
    marginal.push_back(Marginal(0.5 * (data[i] + data[i+1]), //mean
                                1e303, //mode
                                0.5 * (data[i] + data[i+1]), //median
                                std::sqrt((data[i] - data[i+1]) * (data[i] - data[i+1]) /12.), //std dev GUM §4.3.7
                                data[i], //up
                                data[i + 1], //low
                                IdPdf) //same id
                      );
  }
}

Delta::Delta(double val, const std::string &nP, const std::string &id):
  Pdf(nP,id)
{
  std::vector<double> data;
  data.push_back(val);
  setPDF(data);
}

void Delta::addMarginalFromData(const std::vector<double>&data)
{
  for(unsigned int i = 0; i < data.size(); i++)
  {
    marginal.push_back(Marginal(data[i], //mean
                                data[i], //mode
                                data[i], //median
                                0.,      //std dev
                                data[i], //up
                                data[i], //low
                                IdPdf) //same id
                      );
  }
}

LogN::LogN(double val, double f,const std::string &nP,const std::string &id):
  Pdf(nP,id)
{
  std::vector<double> data;
  data.push_back(val);
  data.push_back(f);
  setPDF(data);
}

void LogN::testNumberOfParameters(unsigned int nP)
{
  const std::string method("void LogN::testNumberOfParameters(unsigned int )");
  if(nP % 2 != 0)
        antiochError(method,"I need an even number of parameters for a lognormal pdf (mean and relative uncertainty for all marginals).");
}

void LogN::addMarginalFromData(const std::vector<double>&data)
{

  testNumberOfParameters(data.size());
  if(marginal.empty())
  {
    mu.clear();
    sigma.clear();
  }
  for(unsigned int i = 0; i < data.size(); i += 2)
  {
    mu.push_back(std::log(data[i]));
    sigma.push_back(std::log(data[i+1]));
    marginal.push_back(Marginal(std::exp(mu[i] + sigma[i] * sigma[i]/2.), //mean
                                std::exp(mu[i] - sigma[i] * sigma[i]), //mode
                                std::exp(mu[i]), //median
                                std::sqrt( ( std::exp(sigma[i]*sigma[i]) - 1.)*std::exp(2.*mu[i] + sigma[i]*sigma[i]) ), //std dev
                                data[i] + 4. * std::sqrt( ( std::exp(sigma[i]*sigma[i]) - 1.)*std::exp(2.*mu[i] + sigma[i]*sigma[i]) ),  //up
                                data[i] - 4. * std::sqrt( ( std::exp(sigma[i]*sigma[i]) - 1.)*std::exp(2.*mu[i] + sigma[i]*sigma[i]) ),  //down
                                IdPdf) //same id
                      );
  }
}

LogU::LogU(double up, double down,const std::string &nP,const std::string &id):
  Pdf(nP,id)
{
  std::vector<double> data;
  data.push_back(up);
  data.push_back(down);
  setPDF(data);
}

void LogU::testNumberOfParameters(unsigned int nP)
{
  const std::string method("void LogU::testNumberOfParameters(unsigned int )");
  if(nP % 2 != 0)
        antiochError(method,"I need an even number of parameters for a loguniform pdf (upper and lower bounds for all marginals).");
}

void LogU::addMarginalFromData(const std::vector<double>&data)
{
  testNumberOfParameters(data.size());
  
  if(marginal.empty())
  {
    lnup.clear();
    lndown.clear();
  }
  for(unsigned int i = 0; i < data.size(); i += 2)
  {
    lnup.push_back(data[i]);
    lndown.push_back(data[i+1]);
    marginal.push_back(Marginal((data[i] - data[i+1])/(lnup[i] - lndown[i]), //mean
                                1e303, //mode
                                std::sqrt(data[i] * data[i+1]), //median
                                std::sqrt( ( 
                                        (lnup[i]-lndown[i])*(data[i]*data[i] - data[i+1]*data[i+1]) 
                                        - 2.*(data[i]-data[i+1])*(data[i]-data[i+1]) 
                                      )/
                                      (2.*(lnup[i] - lndown[i])*(lnup[i] - lndown[i]) )  ), 
//std dev VOSE software (http://www.vosesoftware.com/ModelRiskHelp/index.htm#Distributions/Continuous_distributions/LogUniform_distribution.htm)
                                data[i], //up
                                data[i + 1], //low
                                IdPdf) //same id
                      );
  }
}


Diri::Diri(double m, double x, const std::string &nP, const std::string &id):
        Pdf(nP,id),
        x(-1.),gam(-1.),sum(0.)
{
  std::vector<double> data;
  data.push_back(m);
  data.push_back(x);
  setPDF(data);
}

void Diri::addMarginalFromData(const std::vector<double>&data)
{
  const std::string method("void Diri::addMarginalFromData(const std::vector<double>&)");
  if(x == -1.)x = data.back();

  if(marginal.empty())sum = 0.;

  for(unsigned int i = 0; i < data.size() - 1; i++)
  {
    sum += data[i];
    marginal.push_back(Marginal(data[i], //mean
                                -1e303, //mode
                                 1e303, //median
                                -1e303, //std dev
                                1, //up
                                0, //low
                                "beta") //id
                      );
  }
}


double Diri::setBetaMode(double alpha, double beta)
{
  if(alpha > 1. && beta > 1.)
  {
     return ((alpha - 1.)/(alpha + beta -1.));
  }else if(alpha == beta)
  {
     return 0.5;
  }else //no mode, only anti mode, search by hand on the distribution if really needed
  {
     return 1e303;
  }
}

void Diri::initializePdf()
{
  const std::string method("void Diri::initializePdf()");
  if(sum != 1.)antiochError(method,"A Dirichlet pdf must have a fixed sum of one before being initialized.");
  double a(0.),b(0.);
  double gamConstrain(-1.);

  if(x == -1.)//mean relative uncertainty
  {
    x = 0.;
    for(int i = 0; i < getNMarginal(); i++)
    {
       x += getStdDevMar(i)/getMeanMar(i);
    }
    x /= (double)getNMarginal();
  }

  for(int i = 0; i < (int)marginal.size(); i++)
  {
     a += getMeanMar(i) * (1. - getMeanMar(i));
     b += getMeanMar(i) * std::sqrt(a);
     double c;
     (getMeanMar(i) > 0.5)?c = getMeanMar(i):c = 1. - getMeanMar(i);
     if(gamConstrain > c)gamConstrain = c;
  }
  gam = 4./(x*x) * (a/b)*(a/b) - 1.;
  if(gam < gamConstrain)gam = gamConstrain;

  for(int i = 0; i < (int)marginal.size(); i++)
  {
     if(getStdDevMar(i) == -1e303)setStdDevMar(2. * std::sqrt((1./marginal[i].getMean() - 1.)/(gam + 1.)),i);//std from gam (diri) or info (dirg)
     double alpha = 1./getMeanMar(i)*((1./getMeanMar(i) - 1.)/((getMeanMar(i)*getStdDevMar(i))*(getMeanMar(i)*getStdDevMar(i))) -1.);
     double beta = alpha * (1./getMeanMar(i) - 1.);
     setModeMar(setBetaMode(alpha,beta),i);
     setParameterMar(alpha,i); //B(alpha,beta) = B(k,1/theta)
     addParameterMar(beta,i); // B(alpha,beta)
    
/*     setModeMar(setBetaMode(gam * getMeanMar(i),gam * (1. - getMeanMar(i))),i);
     setParameterMar(gam * getMeanMar(i),i); //B(alpha,beta) = B(k,1/theta)
     addParameterMar(gam * (1. - getMeanMar(i)),i); // B(alpha,beta)*/
  }
  
}

DiUn::DiUn(int n, const std::string &nP, const std::string &id):
 Diri(nP,id)
{
  std::vector<double> data;
  double r = 1./(double)n;
  data.insert(data.begin(),n,r);
  data.push_back(1.);
  setPDF(data);
  initializePdf();
}

DirG::DirG(double m, double s,const std::string &nP,const std::string &id):
 Diri(nP,id)
{
   std::vector<double> data;
   data.push_back(m);
   data.push_back(s);
   setPDF(data);
}


void DirG::addMarginalFromData(const std::vector<double>&data)
{
  const std::string method("void DirG::addMarginalFromData(const std::vector<double>&)");

  if(marginal.empty())sum = 0.;

  for(unsigned int i = 0; i < data.size(); i += 2)
  {
    sum += data[i];
    marginal.push_back(Marginal(data[i], //mean
                                -1e303, //mode
                                 1e303, //median
                                data[i+1], //std dev
                                1, //up
                                0, //low
                                "beta") //id
                      );
  }
}


DiUT::DiUT(const std::vector<double> &data,const std::string &nP,const std::string &id):
        DiUn((int)data.size()/2,nP,id)
{
  for(int i = 0; i < getNMarginal(); i++)
  {
     setUpLimitMar(data[2*i],i);
     setDownLimitMar(data[2*i+1],i);
     setMeanMar(1e303,i);
  }
}


DiOr::DiOr(const std::vector<int> &data,const std::string &nP, const std::string &id):
  DiUn((int)data.size(),nP,id)
{
  double n = (double)getNMarginal();
  for(int i = 0; i < getNMarginal(); i++)
  {
     setUpLimitMar(data[2*i],i);
     setMeanMar(1e303,i);
     double r = (double)i;
     setModeMar( (n+1.-r)/(n+(n+1.)/2.) ,i);
  }
}
}
