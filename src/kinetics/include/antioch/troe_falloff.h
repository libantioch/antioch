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

#ifndef ANTIOCH_TROE_FALLOFF_H
#define ANTIOCH_TROE_FALLOFF_H

//Antioch
#include "antioch/math_constants.h"

namespace Antioch
{
  /*!\class TroeFalloff
   *
   * The Troe falloff model is defined by:
   * \f[
   *     \log_{10}\left(F\right) = \log_{10}\left(F_{\text{cent}}\right) / 
   *                               \left[ 
   *                                      1 + \left(
   *                                                 \frac{\log_{10}\left(P_r\right) + c}
   *                                                      {n - d * \left[\log_{10}\left(P_r\right) + c\right]}
   *                                          \right)^2
   *                               \right]
   * \f]
   * with
   * \f[
   *     \begin{array}{r@{\,=\,}l}\toprule
   *     P_r   &   [M] \frac{k_0}{k_\infty}                          \\ 
   *     n     &   0.75 - 1.27 \log_{10}\left(F_{\text{cent}}\right) \\
   *     c     &   0.40 - 0.67 \log_{10}\left(F_{\text{cent}}\right) \\
   *     d     &   0.14                                              \\\bottomrule
   *     \end{array}
   * \f]
   * and
   * \f[
   *     F_{\text{cent}} = (1. - \alpha) * \exp\left(-\frac{T}{T^{***}}\right) + \alpha * \exp\left(-\frac{T}{T^*}\right) + \exp\left(-\frac{T^{**}}{T}\right)
   * \f]
   * 
   * \f$\alpha\f$, \f$T^{*}\f$, \f$T^{**}\f$, \f$T^{***}\f$ being the parameters of the falloff.
   */
  template <typename CoeffType = double>
  class TroeFalloff
  {
  public:
    TroeFalloff(const unsigned int nspec, const CoeffType alpha=0.,
                const CoeffType T3 = 0., const CoeffType T1 = 0.,
                const CoeffType T2 = 1e50);

    ~TroeFalloff();

    template<typename StateType>
    StateType operator()(const StateType& T, const StateType &Pr) const;

    template <typename StateType, typename VectorStateType>
    void F_and_derivatives(const StateType& T, 
                           const StateType &Pr, 
                           const StateType &dPr_dT, 
                           const VectorStateType &dPr_dX,
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
  template<typename StateType>
  inline
  StateType TroeFalloff<CoeffType>::operator()(const StateType& T, const StateType &Pr) const
  {
    StateType Fcent = this->Fcent(T);
    using std::log;

    // c = -0.4 - 0.67 * log10(Fcent)
    // Note log10(x) = (1.0/log(10))*log(x)
    StateType  c = -0.4 - _c_coeff*log(Fcent);

    // n = 0.75 - 1.27 * log10(Fcent)
    // Note log10(x) = (1.0/log(10))*log(x)
    StateType  n = 0.75 - _n_coeff*log(Fcent);

    // Pr = [M] * k0/kinf
    StateType logPr = Constants::log10_to_log<CoeffType>()*log(Pr);

    using std::pow;
    //logF =  log10(Fcent) / [1+((log10(Pr) + c)/(n - 0.14*(log10(Pr) + c) ))^2]
    StateType logF = Constants::log10_to_log<CoeffType>()*log(Fcent)/(1. + pow(((logPr + c)/(n - 0.14*(logPr + c) )),2) );

    return exp(logF);
  }


  template <typename CoeffType>
  template <typename StateType>
  inline
  StateType TroeFalloff<CoeffType>::Fcent(const StateType &T) const
  {
    using std::exp;
     
    StateType one(1.);

    // Fcent = (1.-alpha)*exp(-T/T***) + alpha * exp(-T/T*) + exp(-T**/T)
    StateType Fc = (one - _alpha) * exp(-T/_T3) + _alpha * exp(-T/_T1);

    if(_T2 != 1e50)Fc += exp(-_T2/T);

    return Fc;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void TroeFalloff<CoeffType>::Fcent_and_derivatives(const StateType &T, StateType &Fc, StateType &dFc_dT) const
  {
    using std::exp;
    
    StateType one(1.);

    // Fcent = (1.-alpha)*exp(-T/T***) + alpha * exp(-T/T*) + exp(-T**/T)
    Fc = (one - _alpha) * exp(-T/_T3) + _alpha * exp(-T/_T1);
    dFc_dT = -(one - _alpha)/_T3 * exp(-T/_T3) - _alpha/_T1 * exp(-T/_T1);

    if(_T2 != 1e50)
      {
        Fc += exp(-_T2/T);
        dFc_dT += _T2/(T*T) * exp(-_T2/T);
      }

    return;
  }

  template <typename CoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  void TroeFalloff<CoeffType>::F_and_derivatives(const StateType& T, 
                                                 const StateType &Pr, 
                                                 const StateType &dPr_dT, 
                                                 const VectorStateType &dPr_dX,
                                                 StateType &F,
                                                 StateType &dF_dT,
                                                 VectorStateType &dF_dX) const
  {
    using std::exp;
    using std::log;
    using std::pow;

    antioch_assert_equal_to(dF_dX.size(),this->n_spec);

    //declarations
    VectorStateType dlogPr_dX = Antioch::zero_clone(dF_dX);
    StateType Fcent = Antioch::zero_clone(T);
    StateType dFcent_dT = Antioch::zero_clone(T);
    StateType logFcent = Antioch::zero_clone(T);
    StateType dlogFcent_dT = Antioch::zero_clone(T);
    StateType logPr = Antioch::zero_clone(T);
    StateType dlogPr_dT = Antioch::zero_clone(T);
    StateType upPart = Antioch::zero_clone(T);
    StateType dupPart_dT = Antioch::zero_clone(T);
    StateType downPart = Antioch::zero_clone(T);
    StateType ddownPart_dT = Antioch::zero_clone(T);
    StateType f = Antioch::zero_clone(T);
    StateType df_dT = Antioch::zero_clone(T);
    StateType log10F = Antioch::zero_clone(T);
    StateType logF = Antioch::zero_clone(T);
    StateType dlogF_dT = Antioch::zero_clone(T);

    //params and change of variable
    this->Fcent_and_derivatives(T,Fcent,dFcent_dT);

    logFcent = Constants::log10_to_log<CoeffType>()*log(Fcent);
    dlogFcent_dT = Constants::log10_to_log<CoeffType>()*dFcent_dT/Fcent;
    logPr = Constants::log10_to_log<CoeffType>()*log(Pr);
    dlogPr_dT = Constants::log10_to_log<CoeffType>()*dPr_dT/Pr;

    for(unsigned int ip = 0; ip < dlogPr_dX.size(); ip++)
      {
        dlogPr_dX[ip] = Constants::log10_to_log<CoeffType>()*dPr_dX[ip]/Pr;
      }

    //decomposing the equation
    {
      CoeffType tmp(0.67L);
      upPart = -0.4L - tmp*logFcent + logPr;
      dupPart_dT = -tmp*dlogFcent_dT + dlogPr_dT;
    }

    //dupPart_dX = dlogPr_dX
    {
      CoeffType tmp2(1.1761L);
      CoeffType tmp(0.14L);
      downPart = 0.806 - tmp2*logFcent - tmp*logPr;
      ddownPart_dT = - tmp2*dlogFcent_dT - tmp*dlogPr_dT;
    }

    //ddownPart_dX = -0.14*dlogPr_dX
    f = upPart/downPart;
    df_dT = f * (dupPart_dT/upPart - ddownPart_dT/downPart);
    //df_dX = dlogPr_dX * f * [1./upPart + 0.14/downPart]

    //finally log10F
    {
      CoeffType tmp(2.0L);
      CoeffType tmp2(10.0L);
      CoeffType tmp3(1.0L);
      CoeffType tmp4(0.14L);
      
      log10F = logFcent / (1.L + f*f);
      logF = CoeffType(1.0/std::log10(std::exp(1.0)))*log10F;
      dlogF_dT = dlogFcent_dT / (1.L + f*f) - tmp * df_dT * f * logFcent/( (1.L + f*f) * (1.L + f*f) );

      //dlogF_dX = - 2. * df_dX * f * logFcent/( (1. + f*f) * (1. + f*f) );
      F = exp(logF);
      dF_dT = dlogF_dT * log(tmp2) * F;
      for(unsigned int ip = 0; ip < dF_dX.size(); ip++)
        {
          dF_dX[ip] = - tmp * dlogPr_dX[ip] * f * (tmp3/upPart + tmp4/downPart) * f * logFcent/( (1.L + f*f) * (1.L + f*f) );
        }
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
    _c_coeff( 0.67L*Constants::log10_to_log<CoeffType>() ),
    _n_coeff( 1.27L*Constants::log10_to_log<CoeffType>() )
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
