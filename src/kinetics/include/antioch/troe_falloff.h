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

#ifndef _ANTIOCH_TROE_FALLOFF_H
#define _ANTIOCH_TROE_FALLOFF_H

//Antioch

//C++

namespace Antioch{
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
class TroeFalloff{
     public:
       TroeFalloff(const unsigned int nspec, const CoeffType alpha=0., const CoeffType T3 = 0., const CoeffType T1 = 0., const CoeffType T2 = 1e50);
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

     template <typename StateType>
     StateType _Fcent(const StateType &T) const;
     template <typename StateType>
     void _Fcent_and_derivatives(const StateType &T,
                                StateType &Fc,
                                StateType &dFc_dT) const;
};
  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType TroeFalloff<CoeffType>::operator()(const StateType& T, const StateType &Pr) const
  {
    StateType Fcent = _Fcent(T);
    using std::log10;
// c = -0.4 - 0.67 * log10(Fcent)
    StateType  c = -0.4 - 0.67 * log10(Fcent);
// n = 0.75 - 1.27 * log10(Fcent)
    StateType  n = 0.75 - 1.27 * log10(Fcent);
// Pr = [M] * k0/kinf
    StateType logPr = log10(Pr);

    using std::pow;
//logF =  log10(Fcent) / [1+((log10(Pr) + c)/(n - 0.14*(log10(Pr) + c) ))^2]
    StateType logF = log10(Fcent)/(1. + pow(((logPr + c)/(n - 0.14*(logPr + c) )),2) );

    return exp(logF);
  }


  template <typename CoeffType>
  template <typename StateType>
  inline
  StateType TroeFalloff<CoeffType>::_Fcent(const StateType &T) const
  {
// Fcent = (1.-alpha)*exp(-T/T***) + alpha * exp(-T/T*) + exp(-T**/T)
    using std::exp;
    StateType one(1.);
    StateType Fc = (one - _alpha) * exp(-T/_T3) + _alpha * exp(-T/_T1);
    if(_T2 != 1e50)Fc += exp(-_T2/T);
    return Fc;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void TroeFalloff<CoeffType>::_Fcent_and_derivatives(const StateType &T, StateType &Fc, StateType &dFc_dT) const
  {
// Fcent = (1.-alpha)*exp(-T/T***) + alpha * exp(-T/T*) + exp(-T**/T)
    using std::exp;
    StateType one(1.);
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
    using std::log10;
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
    StateType logF = Antioch::zero_clone(T);
    StateType dlogF_dT = Antioch::zero_clone(T);
//
    StateType tmp = Antioch::zero_clone(T);
    StateType tmp2 = Antioch::zero_clone(T);
    StateType tmp3 = Antioch::zero_clone(T);
    StateType tmp4 = Antioch::zero_clone(T);

//params and change of variable
    _Fcent_and_derivatives(T,Fcent,dFcent_dT);
    logFcent = log10(Fcent);
    tmp = log(10.L);
    dlogFcent_dT = dFcent_dT/(Fcent * tmp);
    logPr = log10(Pr);
    dlogPr_dT = dPr_dT/(Pr * tmp);
    for(unsigned int ip = 0; ip < dlogPr_dX.size(); ip++)dlogPr_dX[ip] = dPr_dX[ip]/(tmp*Pr);    

//decomposing the equation
    tmp = 0.67L;
    upPart = -0.4L - tmp*logFcent + logPr;
    dupPart_dT = -tmp*dlogFcent_dT + dlogPr_dT;
//dupPart_dX = dlogPr_dX
    tmp2 = 1.1761L;
    tmp = 0.14L;
    downPart = 0.806 - tmp2*logFcent - tmp*logPr;
    ddownPart_dT = - tmp2*dlogFcent_dT - tmp*dlogPr_dT;
//ddownPart_dX = -0.14*dlogPr_dX
    f = upPart/downPart;
    df_dT = f * (dupPart_dT/upPart - ddownPart_dT/downPart);
//df_dX = dlogPr_dX * f * [1./upPart + 0.14/downPart]

//finally log10F
    tmp = 2.0L;
    tmp2 = 10.0L;
    tmp3 = 1.0L;
    tmp4 = 0.14L;
    logF = logFcent / (1.L + f*f);
    dlogF_dT = dlogFcent_dT / (1.L + f*f) - tmp * df_dT * f * logFcent/( (1.L + f*f) * (1.L + f*f) );
//dlogF_dX = - 2. * df_dX * f * logFcent/( (1. + f*f) * (1. + f*f) );
    F = pow(tmp,logF);
    dF_dT = dlogF_dT * log(tmp2) * pow(tmp2,logF);
    for(unsigned int ip = 0; ip < dF_dX.size(); ip++)
    {
       dF_dX[ip] = - tmp * dlogPr_dX[ip] * f * (tmp3/upPart + tmp4/downPart) * f * logFcent/( (1.L + f*f) * (1.L + f*f) );
    }
    return;
  }


  template<typename CoeffType>
  inline
  TroeFalloff<CoeffType>::TroeFalloff(const unsigned int nspec, const CoeffType alpha, const CoeffType T3, const CoeffType T1, const CoeffType T2):
    n_spec(nspec),
    _alpha(alpha),
    _T3(T3),
    _T1(T1),
    _T2(T2)
  {
    return;
  }

  template<typename CoeffType>
  inline
  TroeFalloff<CoeffType>::~TroeFalloff()
  {
    return;
  }
}

#endif
