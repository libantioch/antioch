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
 *     \log_{10}\left(F\left) = \log_{10}\left(F_{\text{cent}}\right) / 
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
       TroeFalloff();
       ~TroeFalloff();

     template<typename StateType>
     StateType operator()(const StateType& T, const StateType &Pr) const;

     template <typename StateType, typename VectorStateType>
     StateType F_and_derivatives(const StateType& T, 
                           const StateType &Pr, 
                           const StateType &dPr_dT, 
                           const VectorStateType &dPr_dY,
                           StateType &F,
                           StateType &dF_dT,
                           VectorStateType &dF_dY) const;


     private:
     CoeffType _alpha;
     CoeffType _T3;
     CoeffType _T1;
     CoeffType _T2;

     template <typename StateType>
     StateType _Fcent(const StateType &T);
     template <typename StateType>
     void _Fcent_and_derivatives(const StateType &T,
                                StateType &Fc,
                                StateType &dFc_dT);
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
  StateType TroeFalloff<CoeffType>::_Fcent(const StateType &T)
  {
// Fcent = (1.-alpha)*exp(-T/T***) + alpha * exp(-T/T*) + exp(-T**/T)
    using std::exp;
    StateType Fc = (1. - _alpha) * exp(-T/_T3) + _alpha * exp(-T/_T1);
    if(_T2/T < 1e6)Fc += exp(-_T2/T);
    return Fc;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void TroeFalloff<CoeffType>::_Fcent_and_derivatives(const StateType &T, StateType &Fc, StateType &dFc_dT)
  {
// Fcent = (1.-alpha)*exp(-T/T***) + alpha * exp(-T/T*) + exp(-T**/T)
    using std::exp;
    Fc = (1. - _alpha) * exp(-T/_T3) + _alpha * exp(-T/_T1);
    dFc_dT = -(1. - _alpha)/_T3 * exp(-T/_T3) - _alpha/_T1 * exp(-T/_T1);
    if(_T2/T < 1e6)
    {
      Fc += exp(-_T2/T);
      dFc_dT += _T2/(T*T) * exp(-_T2/T);
    }
    return;
  }

  template <typename CoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  StateType TroeFalloff<CoeffType>::F_and_derivatives(const StateType& T, 
                        const StateType &Pr, 
                        const StateType &dPr_dT, 
                        const VectorStateType &dPr_dY,
                        StateType &F,
                        StateType &dF_dT,
                        VectorStateType &dF_dY) const
  {
    using std::exp;
    using std::log10;
    using std::log;

//params and change of variable
    StateType Fcent,dFcent_dT;
    _Fcent_and_derivatives(T,Fcent,dFcent_dT);
    StateType logFcent = log10(Fcent);
    StateType dlogFcent_dT = dFcent_dT/(Fcent * log(10.));
    StateType logPr = log10(Pr);
    StateType dlogPr_dT = dPr_dT/(Pr * log(10.));
    VectorStateType dlogPr_dY;
    dlogPr_dY.resize(this->n_species(), logPr);
    for(unsigned int ip = 0; ip < dlogPr_dY.size(); ip++)dlogPr_dY[ip] = dPr_dY[ip]/(log(10.)*Pr);    

//decomposing the equation
    StateType upPart = -0.4 - 0.67*logFcent + logPr;
    StateType dupPart_dT = -0.67*dlogFcent_dT + dlogPr_dT;
//dupPart_dY = dlogPr_dY
    StateType downPart = 0.806 - 1.1761*logFcent - 0.14*logPr;
    StateType ddownPart_dT = - 1.1761*dlogFcent_dT - 0.14*dlogPr_dT;
//ddownPart_dY = -0.14*dlogPr_dY
    StateType f = upPart/downPart;
    StateType df_dT = f * (dupPart_dT/upPart - ddownPart_dT/downPart);
//df_dY = dlogPr_dY * f * [1./upPart + 0.14/downPart]

//finally log10F
    StateType logF = logFcent / (1. + f*f);
    StateType dlogF_dT = dlogFcent_dT / (1. + f*f) - 2. * df_dT * f * logFcent/( (1. + f*f) * (1. + f*f) );
//dlogF_dY = - 2. * df_dY * f * logFcent/( (1. + f*f) * (1. + f*f) );
    using std::pow;
    F = pow(10.,logF);
    dF_dT = dlogF_dT * log(10.) * pow(10.,logF);
    dF_dY.resize(this->n_species(), logPr);
    for(unsigned int ip = 0; ip < dF_dY.size(); ip++)
    {
       dF_dY[ip] = - 2. * dlogPr_dY * f * (1./upPart + 0.14/downPart) * f * logFcent/( (1. + f*f) * (1. + f*f) );
    }
    return;
  }


  template<typename CoeffType>
  inline
  TroeFalloff<CoeffType>::TroeFalloff()
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
