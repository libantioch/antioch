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

     CoeffType compute_F(const CoeffType& T, const CoeffType &Pr) const;


     private:
     CoeffType _alpha;
     CoeffType _T3;
     CoeffType _T1;
     CoeffType _T2;
};
  template<typename CoeffType>
  inline
  unsigned int Reaction<CoeffType>::compute_F(const CoeffType& T, const CoeffType &Pr) const
  {
    using std::exp;
// Fcent = (1.-alpha)*exp(-T/T***) + alpha * exp(-T/T*) + exp(-T**/T)
    CoeffType Fcent = (1. - alpha) * exp(-T/T3) + alpha * exp(-T/T1);
    if(T2/T < 1e6)Fcent += exp(-T2/T);

    using std::log10;
// c = -0.4 - 0.67 * log10(Fcent)
    CoeffType  c = -0.4 - 0.67 * log10(Fcent);
// n = 0.75 - 1.27 * log10(Fcent)
    CoeffType  n = 0.75 - 1.27 * log10(Fcent);
// Pr = [M] * k0/kinf
    CoeffType logPr = log10(Pr);

    using std::pow;
//logF = [1+((log10(Pr) + c)/(n - 0.14*(log10(Pr) + c) ))^2 ]^-1 * log10(Fcent)
    CoeffType logF = log10(Fcent)/(1. + pow(((logPr + c)/(n - 0.14*(logPr + c) )),2) );

    return exp(logF);
  }


  template<typename CoeffType>
  inline
  TroeFallOff<CoeffType>::TroeFalloff()
  {
    return;
  }

  template<typename CoeffType>
  inline
  TroeFallOff<CoeffType>::~TroeFalloff()
  {
    return;
  }
}
