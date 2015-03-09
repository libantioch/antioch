//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014 Paul T. Bauman, Benjamin S. Kirk, Sylvain Plessis,
//                    Roy H. Stonger
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

#ifndef ANTIOCH_TROE_FALLOFF_H
#define ANTIOCH_TROE_FALLOFF_H

//Antioch
#include "antioch/math_constants.h"
#include "antioch/cmath_shims.h"

namespace Antioch
{
  /*!\class TroeFalloff
   *
   * The Troe falloff model is defined by:
   * \f[
   *     \log_{10}\left(F\right) = \frac{\log_{10}\left(F_{\text{cent}}\right)}
   *                               {1 + \left[
   *                                          \frac{\log_{10}\left(P_r\right) + c}
   *                                               {n - d \cdot \left[\log_{10}\left(P_r\right) + c\right]}
   *                                     \right]^2}
   * \f]
   * with
   * \f[
   * \begin{split}
   *     P_r   & =  [\mathrm{M}] \frac{k_0}{k_\infty} \\ 
   *     n     & =  0.75 - 1.27 \log_{10}\left(F_{\text{cent}}\right) \\
   *     c     & = - 0.40 - 0.67 \log_{10}\left(F_{\text{cent}}\right) \\
   *     d     & =  0.14                                              \\
   *     F_{\text{cent}} & = (1 - \alpha) \cdot \exp\left(-\frac{T}{T^{***}}\right) + \alpha \cdot \exp\left(-\frac{T}{T^*}\right) + \exp\left(-\frac{T^{**}}{T}\right)
   * \end{split}
   * \f]
   * The derivatives are therefore:
   * \f[
   * \begin{split}
   *     \frac{\partial F_{\text{cent}}}{\partial T} 
   *                           & = \frac{\alpha - 1}{T^{***}} \cdot \exp\left(-\frac{T}{T^{***}}\right) 
   *                            - \frac{\alpha}{T^{*}} \cdot \exp\left(-\frac{T}{T^*}\right) 
   *                            + \frac{T^{**}}{T^{2}} \exp\left(-\frac{T^{**}}{T}\right) \\
   *     \frac{\partial \log_{10}\left(F_\text{cent}\right)}{\partial T} & = \frac{1}{\ln(10) F_\text{cent}}\frac{\partial F_\text{cent}}{\partial T} \\
   *     \frac{\partial n}{\partial T} & = - 1.27 \frac{\partial \log_{10}\left(F_\text{cent}\right)}{\partial T} \\
   *     \frac{\partial c}{\partial T} & = - 0.67 \frac{\partial \log_{10}\left(F_\text{cent}\right)}{\partial T} \\\\
   *     \frac{\partial P_r}{\partial T} & = P_r \left(\frac{\partial k_0}{\partial T} \frac{1}{k_0} - \frac{\partial k_\infty}{\partial T} \frac{1}{k_\infty} \right)\\
   *     \frac{\partial \log_{10}(P_r)}{\partial T} & = \frac{1}{\ln(10) P_r} \frac{\partial P_r}{\partial T} \\\\
   *     \frac{\partial \log_{10}(F)}{\partial T} & = \frac{\partial \log_{10}\left(F_\text{cent}\right)}{\partial T}
   *                                                    \frac{1}{1 + \left[\frac{\log_{10}\left(P_r\right) + c}
   *                                                                            {n - d\left(\log_{10}\left(P_r\right) + c\right)}\right]^2}
   *                                                - \log_{10}\left(F_\text{cent}\right)
   *                                                  2\left[\frac{\log_{10}\left(P_r\right) + c}{n - d \left[\log_{10}\left(P_r\right) + c\right]}\right]^2
   *                                                  \left[\frac{\frac{\partial \log_{10}\left(P_r\right)}{\partial T} + \frac{\partial c}{\partial T}}
   *                                                             {\log_{10}\left(P_r\right) + c}
   *                                                       - \frac{\frac{\partial n}{\partial T} 
   *                                                            - d \left[\frac{\partial \log_{10}\left(P_r\right)}{\partial T} + \frac{\partial c}{\partial T}\right]}
   *                                                              {n - d \left[\log_{10}\left(P_r\right) + c\right]}
   *                                                 \right]
   *                                                          \frac{1}{\left[1 + \left[\frac{\log_{10}\left(P_r\right) + c}
   *                                                                                        {n - d\left(\log_{10}\left(P_r\right) + c\right)}
   *                                                                             \right]^2
   *                                                                  \right]^2}
   *     \\
   *     & = \log_{10}\left(F\right) \left[\frac{\partial \log_{10}\left(F_\text{cent}\right)}{\partial T} \frac{1}{F_\text{cent}}
   *                                      - 2\left[\frac{\log_{10}\left(P_r\right) + c}{n - d \left[\log_{10}\left(P_r\right) + c\right]}\right]^2
   *                                                  \left[\frac{\frac{\partial \log_{10}\left(P_r\right)}{\partial T} + \frac{\partial c}{\partial T}}
   *                                                             {\log_{10}\left(P_r\right) + c}
   *                                                       - \frac{\frac{\partial n}{\partial T} 
   *                                                            - d \left[\frac{\partial \log_{10}\left(P_r\right)}{\partial T} + \frac{\partial c}{\partial T}\right]}
   *                                                              {n - d \left[\log_{10}\left(P_r\right) + c\right]}
   *                                                 \right]
   *                                                          \frac{1}{1 + \left[\frac{\log_{10}\left(P_r\right) + c}
   *                                                                                        {n - d\left(\log_{10}\left(P_r\right) + c\right)}
   *                                                                             \right]^2}
   *       \right] \\
   *     \frac{\partial F}{\partial T} & = \ln(10) F  \frac{\partial \log_{10}\left(F\right)}{\partial T} \\\\\\\\
   *     \frac{\partial P_r}{\partial c_i} & = \frac{k_0}{k_\infty} \\
   *     \frac{\partial \log_{10}(P_r)}{\partial c_i} & = \frac{1}{\ln(10) P_r} \frac{\partial P_r}{\partial c_i} = \frac{1}{\ln(10) [\mathrm{M}]}\\\\
   *     \frac{\partial \log_{10}\left(F\right)}{\partial c_i} & = -\frac{\log_{10}^2\left(F\right)}{\log_{10}\left(F_\text{cent}\right)}
   *                                                               \frac{\partial \log_{10}\left(P_r\right)}{\partial c_i}
   *                                                               \left(1 - \frac{1}{n - d\left[\log_{10}\left(P_r\right) + c\right]}\right)
   *                                                               \left(\log_{10}\left(P_r\right) + c\right) \\
   *     \frac{\partial F}{\partial c_i} & = \ln(10) F \frac{\partial \log_{10}\left(F\right)}{\partial c_i} 
   * \end{split}
   * \f]
   * 
   * \f$\alpha\f$, \f$T^{*}\f$, \f$T^{**}\f$, \f$T^{***}\f$ being the parameters of the falloff.
   */
  template <typename CoeffType = double>
  class TroeFalloff
  {
  public:
    TroeFalloff(const unsigned int nspec, const CoeffType alpha=0,
                const CoeffType T3 = 0, const CoeffType T1 = 0,
                const CoeffType T2 = 1e50);

    ~TroeFalloff();

    void set_alpha(const CoeffType &al);
    void set_T1(const CoeffType &T);
    void set_T2(const CoeffType &T);
    void set_T3(const CoeffType &T);

    template <typename StateType>
    StateType operator()(const StateType &T,
                         const StateType &M,
                         const StateType &k0, 
                         const StateType &kinf) const;

    template <typename StateType, typename VectorStateType>
    void F_and_derivatives(const StateType& T, 
                           const StateType &M,
                           const StateType &k0, 
                           const StateType &dk0_dT, 
                           const StateType &kinf, 
                           const StateType &dkinf_dT, 
                           StateType &F,
                           StateType &dF_dT,
                           VectorStateType &dF_dX) const;

  private:

    unsigned int n_spec;
    CoeffType _alpha;
    CoeffType _T3;
    CoeffType _T1;
    CoeffType _T2;

    //! Precompute coefficient for log conversion.
    /*! This is needed because Eigen doesn't understand log10. */
    CoeffType _c_coeff;

    //! Precompute coefficient for log conversion
    /*! This is needed because Eigen doesn't understand log10. */
    CoeffType _n_coeff;

    template <typename StateType>
    StateType Fcent(const StateType &T) const;

    template <typename StateType>
    void Fcent_and_derivatives( const StateType &T,
                                StateType &Fc,
                                StateType &dFc_dT ) const;

  };

  template<typename CoeffType>
  inline
  void TroeFalloff<CoeffType>::set_alpha(const CoeffType& al)
  {
    _alpha = al;
    return;
  }

  template<typename CoeffType>
  inline
  void TroeFalloff<CoeffType>::set_T1(const CoeffType& T)
  {
    _T1 = T;
    return;
  }

  template<typename CoeffType>
  inline
  void TroeFalloff<CoeffType>::set_T2(const CoeffType& T)
  {
    _T2 = T;
    return;
  }

  template<typename CoeffType>
  inline
  void TroeFalloff<CoeffType>::set_T3(const CoeffType& T)
  {
    _T3 = T;
    return;
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType TroeFalloff<CoeffType>::operator()(const StateType& T,
                                               const StateType &M,
                                               const StateType &k0, 
                                               const StateType &kinf) const
  {

    //compute log(Fcent) once
    StateType logFcent = ant_log(this->Fcent(T));

    // Pr = [M] * k0/kinf
    ANTIOCH_AUTO(StateType) Pr = M * k0/kinf;
    // c = -0.4 - 0.67 * log10(Fcent)
    // Note log10(x) = (1.0/log(10))*log(x)
    StateType  c = - CoeffType(0.4L) - _c_coeff*logFcent;

    // n = 0.75 - 1.27 * log10(Fcent)
    // Note log10(x) = (1.0/log(10))*log(x)
    ANTIOCH_AUTO(StateType) n = CoeffType(0.75L) - _n_coeff*logFcent;
    ANTIOCH_AUTO(StateType) d = constant_clone(T,0.14L);

    StateType log10Pr = Constants::log10_to_log<CoeffType>() * ant_log(Pr);

    //log10F =  log10(Fcent) / [1+((log10(Pr) + c)/(n - d*(log10(Pr) + c) ))^2]
    //logF =  log(Fcent) / [1+((log10(Pr) + c)/(n - d*(log10(Pr) + c) ))^2]
    ANTIOCH_AUTO(StateType) logF =
      logFcent/(1 + ant_pow(((log10Pr + c)/(n - d*(log10Pr + c) )),2) );

    return ant_exp(logF);
  }


  template <typename CoeffType>
  template <typename StateType>
  inline
  StateType TroeFalloff<CoeffType>::Fcent(const StateType &T) const
  {
     
    // Fcent = (1.-alpha)*exp(-T/T***) + alpha * exp(-T/T*) + exp(-T**/T)
    StateType Fc = (1 - _alpha) * ant_exp(-T/_T3) + _alpha * ant_exp(-T/_T1);

    if(_T2 != 1e50)Fc += ant_exp(-_T2/T);

    return Fc;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void TroeFalloff<CoeffType>::Fcent_and_derivatives(const StateType &T, StateType &Fc, StateType &dFc_dT) const
  {
    
    // Fcent = (1.-alpha)*exp(-T/T***) + alpha * exp(-T/T*) + exp(-T**/T)
    Fc = (1 - _alpha) * ant_exp(-T/_T3) + _alpha * ant_exp(-T/_T1);
    dFc_dT = (_alpha - 1)/_T3 * ant_exp(-T/_T3) - _alpha/_T1 * ant_exp(-T/_T1);

    if(_T2 != 1e50)
      {
        Fc += ant_exp(-_T2/T);
        dFc_dT += _T2/(T*T) * ant_exp(-_T2/T);
      }

    return;
  }

  template <typename CoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  void TroeFalloff<CoeffType>::F_and_derivatives(const StateType& T, 
                                                 const StateType &M,
                                                 const StateType &k0, 
                                                 const StateType &dk0_dT, 
                                                 const StateType &kinf, 
                                                 const StateType &dkinf_dT, 
                                                 StateType &F,
                                                 StateType &dF_dT,
                                                 VectorStateType &dF_dX) const
  {

    antioch_assert_equal_to(dF_dX.size(),this->n_spec);

    //declarations
    // Pr and derivatives
    StateType Pr = M * k0/kinf;
    StateType dPr_dT = Pr * (dk0_dT/k0 - dkinf_dT/kinf);
    StateType log10Pr = Constants::log10_to_log<CoeffType>() * ant_log(Pr);
    StateType dlog10Pr_dT = Constants::log10_to_log<CoeffType>()*dPr_dT/Pr;
    VectorStateType dlog10Pr_dX = Antioch::zero_clone(dF_dX);
    for(unsigned int ip = 0; ip < dlog10Pr_dX.size(); ip++)
      {//dlog10Pr_dX = 1/(ln(10)*Pr) * dPr_dX
        dlog10Pr_dX[ip] = Constants::log10_to_log<CoeffType>()/M; //dPr_dX = k0/kinf, Pr = M k0/kinf => dlog10Pr_dX = 1/(ln(10)*M)
      }
    // Fcent and derivatives
    StateType Fcent = Antioch::zero_clone(T);
    StateType dFcent_dT = Antioch::zero_clone(T);
    this->Fcent_and_derivatives(T,Fcent,dFcent_dT);
    StateType dlog10Fcent_dT = Constants::log10_to_log<CoeffType>()*dFcent_dT/Fcent;
    // Compute log(Fcent) once
    StateType logFcent = ant_log(Fcent);
    // n and c and derivatives
    StateType  d = Antioch::constant_clone(T, 0.14L);
    StateType  c = - CoeffType(0.4L) - _c_coeff * logFcent;
    StateType  n = CoeffType(0.75L) - _n_coeff * logFcent;
    StateType dc_dT = - _c_coeff * dFcent_dT/Fcent;
    ANTIOCH_AUTO(StateType) dn_dT = - _n_coeff * dFcent_dT/Fcent;

    //log10F
    StateType logF = logFcent/(1 + ant_pow(((log10Pr + c)/(n - d*(log10Pr + c) )),2));
    StateType dlogF_dT = logF * (dlog10Fcent_dT / Fcent 
                                     - 2 * ant_pow((log10Pr + c)/(n - d * (log10Pr + c)),2)
                                       * ((dlog10Pr_dT + dc_dT)/(log10Pr + c) -
                                          (dn_dT - d * (dlog10Pr_dT + dc_dT))/(n - d * (log10Pr + c))
                                         )
                                       / (1 + ant_pow((log10Pr + c)/(n - d * (log10Pr + c)),2))
                                    );
    VectorStateType dlogF_dX = Antioch::zero_clone(dF_dX);
    for(unsigned int ip = 0; ip < dlog10Pr_dX.size(); ip++)
      {//dlogF_dX = - logF^2/log(Fcent) * dlog10Pr_dX * (1 - 1/(n - d * (log10Pr + c))) * (log10Pr + c)
        dlogF_dX[ip] = - ant_pow(logF,2)/logFcent * dlog10Pr_dX[ip] *(1 - 1/(n - d * (log10Pr + c))) * (log10Pr + c);
      }

    F = ant_exp(logF);
    dF_dT = F * dlogF_dT;
    for(unsigned int ip = 0; ip < dlog10Pr_dX.size(); ip++)
      {//dF_dX =  F * dlogF_dX
        dF_dX[ip] = F * dlogF_dX[ip];
      }

    return;
  }


  template<typename CoeffType>
  inline
  TroeFalloff<CoeffType>::TroeFalloff(const unsigned int nspec, const CoeffType alpha,
                                      const CoeffType T3, const CoeffType T1, const CoeffType T2):
    n_spec(nspec),
    _alpha(alpha),
    _T3(T3),
    _T1(T1),
    _T2(T2),
    _c_coeff( 0.67L * Constants::log10_to_log<CoeffType>() ),
    _n_coeff( 1.27L * Constants::log10_to_log<CoeffType>() )
  {
    return;
  }

  template<typename CoeffType>
  inline
  TroeFalloff<CoeffType>::~TroeFalloff()
  {
    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_TROE_FALLOFF_H
