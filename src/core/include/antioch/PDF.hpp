//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
#ifndef _PDF_DESCRIPTIONS_
#define _PDF_DESCRIPTIONS_

/*!\file Pdf.hpp
 * \brief PDFs definitions
 * \author \me
 *
 * define the pdfs. structure is:
 * a pdf is a std::vector of marginals
 *   - univariate size() = 1
 *   - multivariate size() > 1
 *
 * the marginals contains only the values,
 * the pdf contains the function on
 * how the values relate, the tests
 * and common values
 * (like sum & sum tests for Dirichlet distributions),
 *
 * \todo Documentation on the way...
 */

//Antioch
#include "antioch/Error.hpp"

//C++
#include <vector>
#include <string>

namespace Antioch{

/*!\class marginPdf
 * \brief mother class of all marginals
 *
 * All marginals class MUST derive from it
 * A marginal is a monovariate pdf, the values
 * stored are:
 *   - the mean,
 *   - the median,
 *   - the mode,
 *   - the standard deviation,
 *   - an upper limit,
 *   - a lower limit
 * The id of a marginal is the probability distribution
 * from which the values are derived. In order to be
 * able to store all distribution, a std::vector of double
 * is used to store the possible parameters needed for
 * the distibution (ex. \f$\mathrm{Beta}(\alpha,\beta)\f$).
 *
 * The name is optionnal. It is used to discriminate
 * between the marginals of the same multivariate
 * pdf. Thus giving a name to the marginal enables
 * to tell them easily.
 */
class Marginal
{
  public:
/*!\brief Default constructor
 *
 * It sets all the values to \f$-10^{303}\f$. This
 * value is used to discriminate between initialized
 * and non initialized objects. When the value is not
 * defined, the value \f$10^{303}\f$ should be used.
 */
   Marginal(double mean = -1e303,
            double mode = -1e303,
            double median = -1e303,
            double stdDev = -1e303,
            double upLimit = -1e303,
            double downLimit = -1e303,
            std::string IdMarginal = "No Id",
            std::string nameMarginal = "No name"){}
/*!\brief Default destructor*/
   ~Marginal(){}

/*!\brief Mean setter*/
   void setMean(double nom)                 {mean = nom;}
/*!\brief Mode setter*/
   void setMode(double nom)                 {mode = nom;}
/*!\brief Median setter*/
   void setMedian(double nom)               {median = nom;}
/*!\brief Standard deviation setter*/
   void setStdDev(double unc)               {stdDev = unc;}
/*!\brief Upper limit setter*/
   void setUpLimit(double up)               {upLimit = up;}
/*!\brief Lower limit setter*/
   void setDownLimit(double down)           {downLimit = down;}
/*!\brief Name setter*/
   void setNameMarginal(const std::string& name) {nameMarginal = name;}
/*!\brief ID setter*/
   void setIdMarginal(const std::string& id)     {IdMarginal = id;}

/*!\brief Distribution parameters setter*/
   void setParameterMarginal(const std::vector<double> par) {parameters = par;}
/*!\brief Distribution parameter setter*/
   void setParameterMarginal(double par)               {parameters.clear(); parameters.push_back(par);}
/*!\brief Distribution parameter adder*/
   void addParameterMarginal(double par)               {parameters.push_back(par);}

/*!\brief Mean getter*/
   double getMean()               const {return mean;}
/*!\brief Mode getter*/
   double getMode()               const {return mode;}
/*!\brief Median getter*/
   double getMedian()             const {return median;}
/*!\brief Standard deviation getter*/
   double getStdDev()             const {return stdDev;}
/*!\brief Upper limit getter*/
   double getUpLimit()            const {return upLimit;}
/*!\brief Lower limit getter*/
   double getDownLimit()          const {return downLimit;}
/*!\brief Name getter*/
   const std::string getNameMarginal() const {return nameMarginal;}
/*!\brief ID getter*/
   const std::string getIdMarginal()   const {return IdMarginal;}

/*!\brief Number of parameters of the distribution getter*/
   int getNParametersMarginal()         const {return (int)parameters.size();}
/*!\brief Parameter of the distribution getter*/
   double getParameter(int n)           const;
/*!\brief Parameters of the distribution getter*/
   const std::vector<double> getParameters() const {return parameters;}

  protected:
/*!\brief A double for the mean*/
    double mean;
/*!\brief A double for the mode*/
    double mode;
/*!\brief A double for the median*/
    double median;
/*!\brief A double for the standard Deviation*/
    double stdDev;
/*!\brief A double for the upper limit*/
    double upLimit;
/*!\brief A double for the lower limit*/
    double downLimit;
/*!\brief A std::vector of doubles for the parameters of the distribution*/
    std::vector<double> parameters;
/*!\brief A std::string for the ID*/
    std::string IdMarginal;
/*!\brief A std::string for the name*/
    std::string nameMarginal;
};

/*\class Pdf
 * \brief This class codes for the pdf of the parameters
 *
 * A pdf is a collection of marginals, from one (univariate)
 * to any number (multivariate). This is the mother class,
 * no calculations is performed here, it contains only the
 * setters and getters.
 *
 */
class Pdf
{
  public:
    Pdf(const std::string &nP = "Undefined pdf", const std::string &IDp = "No Id"):namePdf(nP),IdPdf(IDp),nConnected(0){}
    ~Pdf(){}

   Marginal marginalFromIndex(int ind);
   Marginal marginalFromName(const std::string &id) const;
   int getNMarginal() const {return (int)marginal.size();}

   void setMeanMar(double m, int nmar = 0)                {marginal[nmar].setMean(m);}
   void setModeMar(double m, int nmar = 0)                {marginal[nmar].setMode(m);}
   void setMedianMar(double m,int nmar = 0)               {marginal[nmar].setMedian(m);}
   void setStdDevMar(double unc, int nmar = 0)            {marginal[nmar].setStdDev(unc);}
   void setUpLimitMar(double up, int nmar = 0)            {marginal[nmar].setUpLimit(up);}
   void setDownLimitMar(double down, int nmar = 0)        {marginal[nmar].setDownLimit(down);}
   void setNameMar(const std::string& name, int nmar = 0)      {marginal[nmar].setNameMarginal(name);}
   void setIdMar(const std::string& id, int nmar = 0)          {marginal[nmar].setIdMarginal(id);}

   void setParameterMar(const std::vector<double> par,int nmar = 0) {marginal[nmar].setParameterMarginal(par);}
   void setParameterMar(double par, int nmar = 0)              {marginal[nmar].setParameterMarginal(par);}
   void addParameterMar(double par, int nmar = 0)              {marginal[nmar].addParameterMarginal(par);}

   void setMeanMar(double m, const std::string &name)                {marginalFromName(name).setMean(m);}
   void setModeMar(double m, const std::string &name)                {marginalFromName(name).setMode(m);}
   void setMedianMar(double m,const std::string &name)               {marginalFromName(name).setMedian(m);}
   void setStdDevMar(double unc, const std::string &name)            {marginalFromName(name).setStdDev(unc);}
   void setUpLimitMar(double up, const std::string &name)            {marginalFromName(name).setUpLimit(up);}
   void setDownLimitMar(double down, const std::string &name)        {marginalFromName(name).setDownLimit(down);}
   void setIdMarginal(const std::string& name, const std::string &idmar)  {marginalFromName(name).setIdMarginal(idmar);}

   double getMeanMar(int nmar = 0)       const {return marginal[nmar].getMean();}
   double getModeMar(int nmar = 0)       const {return marginal[nmar].getMode();}
   double getMedianMar(int nmar = 0)     const {return marginal[nmar].getMedian();}
   double getStdDevMar(int nmar = 0)     const {return marginal[nmar].getStdDev();}
   double getUpLimitMar(int nmar = 0)    const {return marginal[nmar].getUpLimit();}
   double getDownLimitMar(int nmar = 0)  const {return marginal[nmar].getDownLimit();}
   const std::string getNameMar(int nmar = 0) const {return marginal[nmar].getNameMarginal();}
   const std::string getIdMar(int nmar = 0)   const {return marginal[nmar].getIdMarginal();}

   int getNParametersMar(int nmar = 0)                 const {return marginal[nmar].getNParametersMarginal();}
   double getParameterMar(int n, int nmar = 0)         const {return marginal[nmar].getParameter(n);}
   const std::vector<double> getParametersMar(int nmar = 0) const {return marginal[nmar].getParameters();}

   double getMeanMar(const std::string &name)       const {return marginalFromName(name).getMean();}
   double getModeMar(const std::string &name)       const {return marginalFromName(name).getMode();}
   double getMedianMar(const std::string &name)     const {return marginalFromName(name).getMedian();}
   double getStdDevMar(const std::string &name)     const {return marginalFromName(name).getStdDev();}
   double getUpLimitMar(const std::string &name)    const {return marginalFromName(name).getUpLimit();}
   double getDownLimitMar(const std::string &name)  const {return marginalFromName(name).getDownLimit();}
   const std::string getIdMar(const std::string &name)   const {return marginalFromName(name).getIdMarginal();}

   int getNParametersMar(const std::string &name)                 const {return marginalFromName(name).getNParametersMarginal();}
   double getParameterMar(const std::string &name, int n)         const {return marginalFromName(name).getParameter(n);}
   const std::vector<double> getParametersMar(const std::string &name) const {return marginalFromName(name).getParameters();}

   void setNamePdf(const std::string &name) {namePdf = name;}
   const std::string getNamePdf() const     {return namePdf;}
   void setIdPdf(const std::string &id)     {IdPdf = id;}
   const std::string getIdPdf()   const     {return IdPdf;}

   void addMarginal(const Marginal &member) {marginal.push_back(member);}

   double getMean()                  const {return marginal[0].getMean();}
   const std::vector<double> getMeanVec() const;

   void isItGood();
   bool isItReady(std::string &error);
   void setPDF(const std::vector<double>&data, const std::vector<std::vector<double> > &mat = std::vector<std::vector<double> >());

   void addOneConnection()     {nConnected++;}
   void supressOneConnection() {nConnected--;}
   int nConnection()     const {return nConnected;}

//overload-to-be methods
   virtual void addMarginalFromData(const std::vector<double>&data);
   virtual void initializePdf();

/*
   virtual double getPdf(double x) const;
   virtual double getCdf(double x) const;
*/

  protected:

    virtual void testNumberOfParameters(unsigned int nP);

    std::vector<Marginal> marginal;
    std::vector<std::vector<double> > Cov;
    std::string namePdf;
    std::string IdPdf;

  private:
    int nConnected;
};

/*!\class Norm
 * \brief Normal pdf
 *
 * Three constructor possibilities:
 *   - explicit double for univariate
 *   - std::vector of doubles for univariate
 *   - std::vector of doubles and covariance matrix for multivariate
 *
 * The relations are:
 *   - mode = median = mean
 *   - upLimit   = mean + 4 * sigma
 *   - downLimit = mean - 4 * sigma
 */
class Norm:public Pdf
{
  public:
    Norm(const std::string &nP = "Non initialized Normal pdf", const std::string &id = "norm"):
                Pdf(nP,id){}
    Norm(double m, double stdUnc,const std::string &nP = "Normal pdf", const std::string &id = "norm");
    Norm(const std::vector<double> &data, const std::vector<std::vector<double> > &mat = std::vector<std::vector<double> >(),const std::string &nP = "Normal pdf", const std::string &id = "norm"):
                        Pdf(nP,id)
                        {setPDF(data,mat);}
    ~Norm(){}

    void addMarginalFromData(const std::vector<double>&data);

  protected:
    void testNumberOfParameters(unsigned int nP);
};

class NorT:public Norm
{
  public:
    NorT(const std::string &nP = "Non Initialized Normal truncated pdf", const std::string &id = "nort"):
                Norm(nP,id){}
    NorT(double m, double stdUnc, double up, double down);
    NorT(const std::vector<double> &data, const std::vector<std::vector<double> >&mat = std::vector<std::vector<double> >(),
                const std::string &nP = "Normal truncated pdf", const std::string &id = "nort"):
                Norm(nP,id){setPDF(data,mat);}
    NorT(const std::vector<double> &norm, const std::vector<double> &limits, const std::vector<std::vector<double> >&mat = std::vector<std::vector<double> >(),
                const std::string &nP = "Normal truncated pdf", const std::string &id = "nort");
    ~NorT(){}

    void addMarginalFromData(const std::vector<double>&data);

  protected:
    void testNumberOfParameters(unsigned int nP);
};

class Unif:public Pdf
{
  public:
    Unif(const std::string &nP = "Non initialized Uniform pdf", const std::string &id = "unif"):
                                Pdf(nP,id){}
    Unif(double up,double down,const std::string &nP = "Uniform pdf", const std::string &id = "unif");
    Unif(const std::vector<double> &data, const std::vector<std::vector<double> >&mat = std::vector<std::vector<double> >(),const std::string &nP = "Uniform pdf", const std::string &id = "unif"):
           Pdf(nP,id){setPDF(data,mat);}

    ~Unif(){}

    void addMarginalFromData(const std::vector<double>&data);

  protected:
    void testNumberOfParameters(unsigned int nP);
};

class Delta:public Pdf
{
  public:
    Delta(const std::string &nP = "Non initialized Delta pdf", const std::string &id = "delt"):
                Pdf(nP,id){}
    Delta(double val,const std::string &nP = "Delta pdf", const std::string &id = "delt");
    Delta(const std::vector<double> &val,const std::string &nP = "Delta pdf", const std::string &id = "delt"):
            Pdf(nP,id){setPDF(val);}
    ~Delta(){}

    void addMarginalFromData(const std::vector<double>&data);

};

class LogN:public Pdf
{
  public:
    LogN(const std::string &nP = "Non initialized LogNormale pdf",const std::string &id = "logn"):
                Pdf(nP,id){}
    LogN(double val, double f,const std::string &nP = "LogNormale pdf",const std::string &id = "logn");
    LogN(const std::vector<double> &data, const std::vector<std::vector<double> > &mat = std::vector<std::vector<double> >(),const std::string &nP = "LogNormale pdf",const std::string &id = "logn"):
        Pdf(nP,id){setPDF(data,mat);}
    ~LogN(){}

    void addMarginalFromData(const std::vector<double>&data);

  protected:
    void testNumberOfParameters(unsigned int nP);

  private:
    std::vector<double> mu;
    std::vector<double> sigma;

};

class LogU:public Pdf
{
  public:
    LogU(const std::string &nP = "Non initialized LogUniform pdf",const std::string &id = "logu"):
                        Pdf(nP,id){}
    LogU(double up, double down,const std::string &nP = "LogUniform pdf",const std::string &id = "logu");
    LogU(const std::vector<double> &data, const std::vector<std::vector<double> > &mat = std::vector<std::vector<double> >(),const std::string &nP = "LogUniform pdf",const std::string &id = "logu"):
        Pdf(nP,id){setPDF(data,mat);}
    ~LogU(){}

    void addMarginalFromData(const std::vector<double>&data);

  protected:
    void testNumberOfParameters(unsigned int nP);

  private:
    std::vector<double> lnup;
    std::vector<double> lndown;
};

//from Carrasco07a notations
class Diri:public Pdf
{
  public:
    Diri(const std::string &nP = "Non initialized Dirichlet pdf",const std::string &IDp = "diri"):
        Pdf(nP,IDp),
        x(-1.),gam(-1.),sum(0.){}
    Diri(double m, double x, const std::string &nP = "Dirichlet pdf", const std::string &id = "diri");
    Diri(const std::vector<double> &data, const std::string &nP = "Dirichlet pdf", const std::string &id = "diri"):
        Pdf(nP,id),
        x(-1.),gam(-1.),sum(0.)
                {setPDF(data);initializePdf();}
    ~Diri(){}

    virtual void initializePdf();

    virtual void addMarginalFromData(const std::vector<double>&data);

  protected:
    double setBetaMode(double alpha, double beta);

    double x;
    double gam;
    double sum;

};

class DiUn:public Diri //DiUn = Diri(1,1,...,1)
{
  public:
    DiUn(const std::string &nP = "Non initialized Dirichlet Uniform pdf",const std::string &id = "diun"):
                Diri(nP,id){}
    DiUn(int n,const std::string &nP = "Dirichlet Uniform pdf",const std::string &id = "diun");
    ~DiUn(){}

};

class DirG:public Diri
{
  public:
    DirG(const std::string &nP = "Non initialized Dirichlet Generilized pdf",const std::string &id = "dirg"):
                Diri(nP,id){}
    DirG(double m, double stdUnc,const std::string &nP = "Dirichlet Generilized pdf",const std::string &id = "dirg");
    DirG(const std::vector<double> &data,const std::string &nP = "Dirichlet Generilized pdf",const std::string &id = "dirg"):
        Diri(nP,id){setPDF(data);initializePdf();}
    ~DirG(){}

    void addMarginalFromData(const std::vector<double>&data);

};

class DiUT:public DiUn
{
  public:
    DiUT(const std::string &nP = "Non initialized Dirichlet Uniform Troncated pdf",const std::string &id = "diut"):
                DiUn(nP,id){}
    DiUT(double up, double down,const std::string &nP = "Dirichlet Uniform Troncated pdf",const std::string &id = "diut"):
        DiUn(nP,id)
           {antiochError("DiUT::DiUT(double)","Please don't use this constructor, use DiUT::DiUT(const std::vector<double>&)");}
    DiUT(const std::vector<double> &data,const std::string &nP = "Dirichlet Uniform Troncated pdf",const std::string &id = "diut");
    ~DiUT(){}

};

class DiOr:public DiUn
{
  public:
    DiOr(const std::string &nP = "Non initialized Dirichlet Uniform Orderer pdf", const std::string &id = "dior"):
                DiUn(nP,id){}
    DiOr(double ind,const std::string &nP = "Dirichlet Uniform Orderer pdf", const std::string &id = "dior"):
                DiUn(nP,id)
                {antiochError("DiOr::DiOr(double)","Please don't use this constructor, use DiOr::DiOr(const std::vector<int>&)");}
    DiOr(const std::vector<int> &data,const std::string &nP = "Dirichlet Uniform Orderer pdf", const std::string &id = "dior");
    ~DiOr(){}

};

/*!\brief Id of the pdfs*/
const std::string KNOWN_PDF[] = {"norm","nort","delt","unif","logn","logu","diri","dirg","diut","diun","dior"};
/*!\brief Number of known pdfs*/
const int nPDF = sizeof(KNOWN_PDF)/sizeof(KNOWN_PDF[0]);
/*!\brief Fork function to obtain the correct pdf object pointer*/
Pdf *objectPDF(const std::string &type,
               const std::vector<double> &parameters = std::vector<double>(), 
               const std::vector<std::vector<double> > &mat = std::vector<std::vector<double> >());

}
#endif
