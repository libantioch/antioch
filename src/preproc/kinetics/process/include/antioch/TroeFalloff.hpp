//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
#ifndef _TROE_FALLOFF_
#define _TROE_FALLOFF_

#include "antioch/FalloffProcess.hpp"

namespace Antioch{
/*!\file TroeFalloff.hpp
 * \brief Contains the TroeFalloff class
 *
 * \class TroeFalloff
 * \brief Troe form of falloff
 *
 * Troe expression for the \f$F\f$ expression
 * is given by the relation:
 * \f[
 *   \log\left(F\right) = \left[1 + \left(\frac{\log\left(P_r\right) + c}
 *                                             {n - d\left[\log\left(P_r\right) + c\right]}\right)^2
 *                       \right]^{-1} \log\left(F_{\text{cent}}\right)
 * \f]
 * with 
 * \f[
 *    \begin{array}{r@{\;=\;}l}
 *    c & -0.4 - 0.67 \log\left(F_{\text{cent}}\right) \\
 *    n & 0.75 - 1.27\log\left(F_{\text{cent}}\right) \\
 *    d & 0.14 \\
 *    \end{array}
 * \f]
 * and \f$F_{\text{cent}}\f$ defined by
 * \f[
 *   F_\text{cent} =
 *    (1-\alpha)\exp\left(-\frac{T}{T^{***}}\right) + \alpha \exp\left(-\frac{T}{T^*}\right) + \exp\left(-\frac{T^{**}}{T}\right)
 * \f]
 * \f$\alpha\f$, \f$T^*\f$, \f$T^{**}\f$ and \f$T^{***}\f$ being parameters,
 * respectively with no unit and K.
 */
class TroeFalloff:public FalloffProcess
{
   public:
/*!\brief Default constructor*/
      TroeFalloff():T(NULL),
                    alpha("alpha","",CORE_UNCERTAINTY_TYPE_NONE),
                    T1("T*","K",CORE_UNCERTAINTY_TYPE_NONE),
                    T2("T**","K",CORE_UNCERTAINTY_TYPE_NONE),
                    T3("T***","K",CORE_UNCERTAINTY_TYPE_NONE),
                    alphaCal("alpha","",CORE_UNCERTAINTY_TYPE_NONE),
                    T1Cal("T*","K",CORE_UNCERTAINTY_TYPE_NONE),
                    T2Cal("T**","K",CORE_UNCERTAINTY_TYPE_NONE),
                    T3Cal("T***","K",CORE_UNCERTAINTY_TYPE_NONE){setKineticsProcess("Troe falloff");}
/*!\brief Copy constructor*/
      TroeFalloff(const TroeFalloff &rhs);
/*!\brief Constructor*/
      TroeFalloff(const std::string &kinMod, 
                         const std::vector<ParameterPhy> &pars0, 
                         const std::vector<ParameterPhy> &parsinf, 
                         const ParameterPhy &al  = ParameterPhy("alpha","",CORE_UNCERTAINTY_TYPE_NONE),
                         const ParameterPhy &To  = ParameterPhy("T*"   ,"K",CORE_UNCERTAINTY_TYPE_NONE),
                         const ParameterPhy &Tt  = ParameterPhy("T**"  ,"K",CORE_UNCERTAINTY_TYPE_NONE),
                         const ParameterPhy &Tth = ParameterPhy("T***" ,"K",CORE_UNCERTAINTY_TYPE_NONE),
                         ParameterPhy *m = NULL);
/*!\brief Destructor*/
      ~TroeFalloff(){}

     void setAlpha(const ParameterPhy &al);
     void setT1(const ParameterPhy &t);
     void setT2(const ParameterPhy &t);
     void setT3(const ParameterPhy &t);

     void setAlpha(double al) {setAlpha(ParameterPhy("alpha",al,0.,CORE_UNCERTAINTY_TYPE_NONE,""));}
     void setT1(double t)     {setT1(ParameterPhy("T*",t,0.,CORE_UNCERTAINTY_TYPE_NONE,"K"));}
     void setT2(double t)     {setT2(ParameterPhy("T**",t,0.,CORE_UNCERTAINTY_TYPE_NONE,"K"));}
     void setT3(double t)     {setT3(ParameterPhy("T***",t,0.,CORE_UNCERTAINTY_TYPE_NONE,"K"));}

/*
 * names are mandatorily 
 * alpha
 * T*
 * T**
 * T***
 */
     void setFalloffParameters(const std::vector<ParameterPhy> &pars);
/*!\brief All parameters getter*/
     const std::vector<ParameterPhy> getParameters() const;

     const ParameterPhy getAlpha() const {return alpha;}
     const ParameterPhy getT1()    const {return T1;}
     const ParameterPhy getT2()    const {return T2;}
     const ParameterPhy getT3()    const {return T3;}

/*!\brief Default configuration alpha: unit*/
     const std::string getDefaultAlphaUnit()       const {return std::string();}
/*!\brief Default configuration alpha: min*/
     double getDefaultAlphaMin()                   const {return 0.;}
/*!\brief Default configuration alpha: max*/
     double getDefaultAlphaMax()                   const {return 1.;}
/*!\brief Default configuration alpha: all in ParameterPhy*/
     const ParameterPhy getDefaultAlphaParameter() const;
/*!\brief Default configuration T*: unit*/
     const std::string getDefaultT1Unit()          const {return "K";}
/*!\brief Default configuration T*: min*/
     double getDefaultT1Min()                      const {return 0.;}
/*!\brief Default configuration T*: max*/
     double getDefaultT1Max()                      const {return 2000.;}
/*!\brief Default configuration T*: all in ParameterPhy*/
     const ParameterPhy getDefaultT1Parameter()    const;
/*!\brief Default configuration T**: unit*/
     const std::string getDefaultT2Unit()          const {return "K";}
/*!\brief Default configuration T**: min*/
     double getDefaultT2Min()                      const {return 0.;}
/*!\brief Default configuration T**: max*/
     double getDefaultT2Max()                      const {return 1e6;}
/*!\brief Default configuration T**: all in ParameterPhy*/
     const ParameterPhy getDefaultT2Parameter()    const;
/*!\brief Default configuration T***: unit*/
     const std::string getDefaultT3Unit()          const {return "K";}
/*!\brief Default configuration T***: min*/
     double getDefaultT3Min()                      const {return 0.;}
/*!\brief Default configuration T***: max*/
     double getDefaultT3Max()                      const {return 2000.;}
/*!\brief Default configuration T***: all in ParameterPhy*/
     const ParameterPhy getDefaultT3Parameter()    const;



/*!\brief Default configuration: unit*/
     const std::string getDefaultFUnit(int i)            const;
/*!\brief Default configuration: min*/
     double getDefaultFMin(int i)                        const;
/*!\brief Default configuration: max*/
     double getDefaultFMax(int i)                        const;
/*!\brief Default configuration: all in ParameterPhy*/
     const ParameterPhy getDefaultFPriorParameter(int i) const;

/*!\brief Default prior configuration: min value*/
     double getDefaultPriorMinValue(int i)              const;
/*!\brief Default prior configuration: max value*/
     double getDefaultPriorMaxValue(int i)              const;
/*!\brief Default prior configuration: unit*/
     const std::string getDefaultPriorUnit(int i)       const;
/*!\brief Default prior configuration: parameter*/
     const ParameterPhy getDefaultPriorParameter(int i) const;


/*!\brief Number of parameters needed for falloff (F)*/
     unsigned int nFPars()        const {return 3;}
/*!\brief Number of optionnal parameters needed for falloff (F)*/
     unsigned int nFOptionPars()  const {return 1;} //T2 sometimes ignored

     void setTemperature(ParameterPhy *Tptr);

   private:

     /*\brief reset the parameters for F*/
     void resetF(std::vector<double> &pars);
     /*\brief reset Alpha*/
     void resetAlpha(double al)                         {resetPar(alpha,al);}
     /*\brief reset T1*/
     void resetT1(double t)                             {resetPar(T1,t);}
     /*\brief reset T2*/
     void resetT2(double t)                             {resetPar(T2,t);}
     /*\brief reset T3*/
     void resetT3(double t)                             {resetPar(T3,t);}
//F
     ParameterPhy FT(double t, int i = 0)               {return exp(logFT(t,i));}
     ParameterPhy FM(double m, int i = 0)               {return exp(logFM(m,i));}
     double FTM(double t, double m, int i = 0)    const {return std::exp(logFTM(t,m,i));}

//logF
     ParameterPhy logFT(double t, int i = 0);
     ParameterPhy logFM(double m, int i = 0);
     double logFTM(double t, double m, int i = 0) const;

//c
     ParameterPhy c(int i = 0)                          {return CEquation<ParameterPhy>(logFcent(i));}
     double cT(double t, int i = 0)               const {return CEquation<double>(logFcentT(t,i));}

//n
     ParameterPhy n(int i = 0)                          {return NEquation<ParameterPhy>(logFcent(i));}
     double nT(double t, int i = 0)               const {return NEquation<double>(logFcentT(t,i));}

//Fcent
     ParameterPhy logFcent(int nF = 0);
     double logFcentT(double t, int i = 0)        const;

//templated equations
     template <class P>
     P logFEquation(const P& Fc, const P& pr, const P& C, const P& N) const;
     template <class P>
     P logFcentEquation(const P& temp, const P& t3, const P& t1, const P& al, const P& t2) const;
     template <class P>
     P CEquation(const P& Fc) const;
     template <class P>
     P NEquation(const P& Fc) const;

     ParameterPhy *T;
     static double d;

     ParameterPhy alpha;
     ParameterPhy T1;
     ParameterPhy T2;
     ParameterPhy T3;

     ParameterPhy alphaCal;
     ParameterPhy T1Cal;
     ParameterPhy T2Cal;
     ParameterPhy T3Cal;
};
}
#endif
