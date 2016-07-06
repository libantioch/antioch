//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014-2016 Paul T. Bauman, Benjamin S. Kirk,
//                         Sylvain Plessis, Roy H. Stonger
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
//--------------------------------------------------------------------------

#ifndef ANTIOCH_STOCKMAYER_POTENTIAL_H
#define ANTIOCH_STOCKMAYER_POTENTIAL_H

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/metaprogramming_decl.h"
#include "antioch/cmath_shims.h"

//C++

namespace Antioch{


  /* Stockmayer potential from
     Monchich & Mason 1961, surface[T][delta]

     In case of temperatures that are too high
     (\f$T > 100 \eps_{\mathrm{LJ}}\f$), we extrapolate using the
     following expression:
     \f[
         \left<\Omega^{(l,s)*}\right> \proto \sqrt[6]{T^{*}}
     \f]
     The proportionality factor is found using the last data point at \f$T^* = 100\f$

     Two points are added when extrapolating to \f$T^* + \Delta T^*\f$ is required,
     \f$100 + \Delta T^*\f$ and \f$100 + frac{\Delta T^*}{2}\f$.
   */
  template <typename CoeffType> //typename VectorCoeffType, typename MatrixCoeffType>
  class StockmayerPotential
  {
        public:
          StockmayerPotential();
          //! In case of temperature extrapolation, T is tested
          template<typename StateType>
          StockmayerPotential(const StateType & T);
          ~StockmayerPotential();

          const std::vector<CoeffType> temperature()     const {return _T;}

          const std::vector<CoeffType> log_temperature() const {return _logT;}

          const std::vector<CoeffType> delta()           const {return _delta;}

          const std::vector<std::vector<CoeffType> > omega_1_1() const {return _omega_1_1;}
          const std::vector<std::vector<CoeffType> > omega_2_2() const {return _omega_2_2;}

          /*! Extrapolation*/
          template <typename StateType>
          void extrapolate_to(const StateType &T);

          CoeffType max_reduced_temperature() const {return _max_reduced_T;}

        private:

          void init();


          template <typename Scalar>
          const Scalar omega_2_2(const Scalar & T_red, const Scalar & dipole_red) const {antioch_not_implemented();}

          template <typename Scalar>
          const Scalar omega_1_1(const Scalar & T_red, const Scalar & dipole_red) const {antioch_not_implemented();}

                unsigned int _T_size;
          const unsigned int _delta_size; // never changes, no extrapolation (yet?) in this direction
//
          std::vector<CoeffType>               _T, _logT;
          std::vector<CoeffType>               _delta;
          std::vector<std::vector<CoeffType> > _omega_1_1;
          std::vector<std::vector<CoeffType> > _omega_2_2;

          const CoeffType _max_reduced_T;

  };

  template <typename CoeffType>
  template <typename StateType>
  inline
  StockmayerPotential<CoeffType>::StockmayerPotential(const StateType & T):
        _T_size(37),
        _delta_size(8),
        _T(_T_size,0.),_logT(_T_size,0.),
        _delta(_delta_size,0.),
        _omega_1_1(_T_size,std::vector<CoeffType>(_delta_size,0.)),
        _omega_2_2(_T_size,std::vector<CoeffType>(_delta_size,0.)),
        _max_reduced_T(100.)
  {
     this->init();

      // tested here
     this->extrapolate_to(T);

     return;
  }

  template <typename CoeffType>
  inline
  StockmayerPotential<CoeffType>::StockmayerPotential():
        _T_size(37),
        _delta_size(8),
        _T(_T_size,0.),_logT(_T_size,0.),
        _delta(_delta_size,0.),
        _omega_1_1(_T_size,std::vector<CoeffType>(_delta_size,0.)),
        _omega_2_2(_T_size,std::vector<CoeffType>(_delta_size,0.)),
        _max_reduced_T(100.)
  {
     this->init();

     return;
  }

  template <typename CoeffType>
  inline
  StockmayerPotential<CoeffType>::~StockmayerPotential()
  {
     return;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void StockmayerPotential<CoeffType>::extrapolate_to(const StateType & T)
  {
      // !\todo better vectorization management
      //
      //  we need to have something as a raw_max,
      //  probably in the form of rawify<StateType,Operation>()
      //  to avoid repetition. Note that two non equivalent
      //  strategies are possible here:
      //     - VectorType<raw_type<StateType> >: all raw_values are equivalent, looking for global Operation
      //     - Operation<StateType>() at each level, looking for local Operation at every level
      CoeffType Delta_T = (max(T) - _max_reduced_T);

        // extrapolation required
      if(Delta_T > 0)
      {

        //resize
        _omega_1_1.push_back(std::vector<CoeffType>(_delta_size)); // + Delta T/2
        _omega_1_1.push_back(std::vector<CoeffType>(_delta_size)); // + Delta T
        _omega_2_2.push_back(std::vector<CoeffType>(_delta_size)); // + Delta T/2
        _omega_2_2.push_back(std::vector<CoeffType>(_delta_size)); // + Delta T

        _T_size += 2;
        _T.push_back(_max_reduced_T + Delta_T/2.);
        _T.push_back(_max_reduced_T + Delta_T);
        _logT.push_back(ant_log(_T[_T_size - 2]));
        _logT.push_back(ant_log(_T[_T_size - 1]));

        for(unsigned int d = 0; d < _delta_size; d++)
        {
           CoeffType K = _omega_1_1[_T_size - 3][d] / (ant_pow(_T[_T_size - 3],CoeffType(1.L/6.L)));
           _omega_1_1[_T_size - 2][d] = K * _T[_T_size - 2];
           _omega_1_1[_T_size - 1][d] = K * _T[_T_size - 1];

           K = _omega_2_2[_T_size - 3][d] / (ant_pow(_T[_T_size - 3],CoeffType(1.L/6.L)));
           _omega_2_2[_T_size - 2][d] = K * _T[_T_size - 2];
           _omega_2_2[_T_size - 1][d] = K * _T[_T_size - 1];
        }
     }
  }

  template <typename CoeffType>
  inline
  void StockmayerPotential<CoeffType>::init()
  {
     ///delta
     _delta[0] = 0.;
     _delta[1] = 0.25;
     _delta[2] = 0.5;
     _delta[3] = 0.75;
     _delta[4] = 1.0;
     _delta[5] = 1.5;
     _delta[6] = 2.0;
     _delta[7] = 2.5;

    ///T
     _T[0]  = (0.1);    _logT[0]  = ant_log(_T[0]);
     _T[1]  = (0.2);    _logT[1]  = ant_log(_T[1]);
     _T[2]  = (0.3);    _logT[2]  = ant_log(_T[2]);
     _T[3]  = (0.4);    _logT[3]  = ant_log(_T[3]);
     _T[4]  = (0.5);    _logT[4]  = ant_log(_T[4]);
     _T[5]  = (0.6);    _logT[5]  = ant_log(_T[5]);
     _T[6]  = (0.7);    _logT[6]  = ant_log(_T[6]);
     _T[7]  = (0.8);    _logT[7]  = ant_log(_T[7]);
     _T[8]  = (0.9);    _logT[8]  = ant_log(_T[8]);
     _T[9]  = (1.0);    _logT[9]  = ant_log(_T[9]);
     _T[10] = (1.2);    _logT[10] = ant_log(_T[10]);
     _T[11] = (1.4);    _logT[11] = ant_log(_T[11]);
     _T[12] = (1.6);    _logT[12] = ant_log(_T[12]);
     _T[13] = (1.8);    _logT[13] = ant_log(_T[13]);
     _T[14] = (2.0);    _logT[14] = ant_log(_T[14]);
     _T[15] = (2.5);    _logT[15] = ant_log(_T[15]);
     _T[16] = (3.0);    _logT[16] = ant_log(_T[16]);
     _T[17] = (3.5);    _logT[17] = ant_log(_T[17]);
     _T[18] = (4.0);    _logT[18] = ant_log(_T[18]);
     _T[19] = (5.0);    _logT[19] = ant_log(_T[19]);
     _T[20] = (6.0);    _logT[20] = ant_log(_T[20]);
     _T[21] = (7.0);    _logT[21] = ant_log(_T[21]);
     _T[22] = (8.0);    _logT[22] = ant_log(_T[22]);
     _T[23] = (9.0);    _logT[23] = ant_log(_T[23]);
     _T[24] = (10.0);   _logT[24] = ant_log(_T[24]);
     _T[25] = (12.0);   _logT[25] = ant_log(_T[25]);
     _T[26] = (14.0);   _logT[26] = ant_log(_T[26]);
     _T[27] = (16.0);   _logT[27] = ant_log(_T[27]);
     _T[28] = (18.0);   _logT[28] = ant_log(_T[28]);
     _T[29] = (20.0);   _logT[29] = ant_log(_T[29]);
     _T[30] = (25.0);   _logT[30] = ant_log(_T[30]);
     _T[31] = (30.0);   _logT[31] = ant_log(_T[31]);
     _T[32] = (35.0);   _logT[32] = ant_log(_T[32]);
     _T[33] = (40.0);   _logT[33] = ant_log(_T[33]);
     _T[34] = (50.0);   _logT[34] = ant_log(_T[34]);
     _T[35] = (75.0);   _logT[35] = ant_log(_T[35]);
     _T[36] = (100.0);  _logT[36] = ant_log(_T[36]);

     //// <omega(1,1)*>
     _omega_1_1[0][0]  = 4.0079;  _omega_1_1[0][1]  = 4.002;  _omega_1_1[0][2]  = 4.655;  _omega_1_1[0][3]  = 4.521;  _omega_1_1[0][4]  = 6.454;  _omega_1_1[0][5]  = 8.214;  _omega_1_1[0][6]  = 9.824;  _omega_1_1[0][7]  = 11.31; 
     _omega_1_1[1][0]  = 3.1300;  _omega_1_1[1][1]  = 3.164;  _omega_1_1[1][2]  = 3.355;  _omega_1_1[1][3]  = 3.721;  _omega_1_1[1][4]  = 4.198;  _omega_1_1[1][5]  = 5.230;  _omega_1_1[1][6]  = 6.225;  _omega_1_1[1][7]  = 7.160; 
     _omega_1_1[2][0]  = 2.6494;  _omega_1_1[2][1]  = 2.657;  _omega_1_1[2][2]  = 2.770;  _omega_1_1[2][3]  = 3.002;  _omega_1_1[2][4]  = 3.319;  _omega_1_1[2][5]  = 4.054;  _omega_1_1[2][6]  = 4.785;  _omega_1_1[2][7]  = 5.483; 
     _omega_1_1[3][0]  = 2.3144;  _omega_1_1[3][1]  = 2.320;  _omega_1_1[3][2]  = 2.402;  _omega_1_1[3][3]  = 2.572;  _omega_1_1[3][4]  = 2.812;  _omega_1_1[3][5]  = 3.386;  _omega_1_1[3][6]  = 3.972;  _omega_1_1[3][7]  = 4.539; 
     _omega_1_1[4][0]  = 2.0661;  _omega_1_1[4][1]  = 2.073;  _omega_1_1[4][2]  = 2.140;  _omega_1_1[4][3]  = 2.278;  _omega_1_1[4][4]  = 2.472;  _omega_1_1[4][5]  = 2.946;  _omega_1_1[4][6]  = 3.437;  _omega_1_1[4][7]  = 3.918; 
     _omega_1_1[5][0]  = 1.8767;  _omega_1_1[5][1]  = 1.885;  _omega_1_1[5][2]  = 1.944;  _omega_1_1[5][3]  = 2.060;  _omega_1_1[5][4]  = 2.225;  _omega_1_1[5][5]  = 2.628;  _omega_1_1[5][6]  = 3.054;  _omega_1_1[5][7]  = 3.474; 
     _omega_1_1[6][0]  = 1.7293;  _omega_1_1[6][1]  = 1.738;  _omega_1_1[6][2]  = 1.791;  _omega_1_1[6][3]  = 1.893;  _omega_1_1[6][4]  = 2.036;  _omega_1_1[6][5]  = 2.388;  _omega_1_1[6][6]  = 2.273;  _omega_1_1[6][7]  = 3.137; 
     _omega_1_1[7][0]  = 1.6122;  _omega_1_1[7][1]  = 1.622;  _omega_1_1[7][2]  = 1.670;  _omega_1_1[7][3]  = 1.760;  _omega_1_1[7][4]  = 1.886;  _omega_1_1[7][5]  = 2.198;  _omega_1_1[7][6]  = 2.535;  _omega_1_1[7][7]  = 2.872; 
     _omega_1_1[8][0]  = 1.5175;  _omega_1_1[8][1]  = 1.527;  _omega_1_1[8][2]  = 1.572;  _omega_1_1[8][3]  = 1.653;  _omega_1_1[8][4]  = 1.765;  _omega_1_1[8][5]  = 2.044;  _omega_1_1[8][6]  = 2.349;  _omega_1_1[8][7]  = 2.657; 
     _omega_1_1[9][0]  = 1.4398;  _omega_1_1[9][1]  = 1.450;  _omega_1_1[9][2]  = 1.490;  _omega_1_1[9][3]  = 1.564;  _omega_1_1[9][4]  = 1.665;  _omega_1_1[9][5]  = 1.917;  _omega_1_1[9][6]  = 2.196;  _omega_1_1[9][7]  = 2.478; 
     _omega_1_1[10][0] = 1.3204;  _omega_1_1[10][1] = 1.330;  _omega_1_1[10][2] = 1.364;  _omega_1_1[10][3] = 1.425;  _omega_1_1[10][4] = 1.509;  _omega_1_1[10][5] = 1.720;  _omega_1_1[10][6] = 1.956;  _omega_1_1[10][7] = 2.199; 
     _omega_1_1[11][0] = 1.2336;  _omega_1_1[11][1] = 1.242;  _omega_1_1[11][2] = 1.272;  _omega_1_1[11][3] = 1.324;  _omega_1_1[11][4] = 1.394;  _omega_1_1[11][5] = 1.573;  _omega_1_1[11][6] = 1.777;  _omega_1_1[11][7] = 1.990; 
     _omega_1_1[12][0] = 1.1679;  _omega_1_1[12][1] = 1.176;  _omega_1_1[12][2] = 1.202;  _omega_1_1[12][3] = 1.246;  _omega_1_1[12][4] = 1.306;  _omega_1_1[12][5] = 1.461;  _omega_1_1[12][6] = 1.639;  _omega_1_1[12][7] = 1.827; 
     _omega_1_1[13][0] = 1.1166;  _omega_1_1[13][1] = 1.124;  _omega_1_1[13][2] = 1.146;  _omega_1_1[13][3] = 1.185;  _omega_1_1[13][4] = 1.237;  _omega_1_1[13][5] = 1.372;  _omega_1_1[13][6] = 1.530;  _omega_1_1[13][7] = 1.698; 
     _omega_1_1[14][0] = 1.0753;  _omega_1_1[14][1] = 1.082;  _omega_1_1[14][2] = 1.102;  _omega_1_1[14][3] = 1.135;  _omega_1_1[14][4] = 1.181;  _omega_1_1[14][5] = 1.300;  _omega_1_1[14][6] = 1.441;  _omega_1_1[14][7] = 1.592; 
     _omega_1_1[15][0] = 1.0006;  _omega_1_1[15][1] = 1.005;  _omega_1_1[15][2] = 1.020;  _omega_1_1[15][3] = 1.046;  _omega_1_1[15][4] = 1.080;  _omega_1_1[15][5] = 1.170;  _omega_1_1[15][6] = 1.278;  _omega_1_1[15][7] = 1.397; 
     _omega_1_1[16][0] = 0.95003; _omega_1_1[16][1] = 0.9538; _omega_1_1[16][2] = 0.9656; _omega_1_1[16][3] = 0.9852; _omega_1_1[16][4] = 1.012;  _omega_1_1[16][5] = 1.082;  _omega_1_1[16][6] = 1.168;  _omega_1_1[16][7] = 1.265; 
     _omega_1_1[17][0] = 0.91311; _omega_1_1[17][1] = 0.9162; _omega_1_1[17][2] = 0.9256; _omega_1_1[17][3] = 0.9413; _omega_1_1[17][4] = 0.9626; _omega_1_1[17][5] = 1.019;  _omega_1_1[17][6] = 1.090;  _omega_1_1[17][7] = 1.170; 
     _omega_1_1[18][0] = 0.88453; _omega_1_1[18][1] = 0.8871; _omega_1_1[18][2] = 0.8948; _omega_1_1[18][3] = 0.9076; _omega_1_1[18][4] = 0.9252; _omega_1_1[18][5] = 0.9721; _omega_1_1[18][6] = 1.031;  _omega_1_1[18][7] = 1.098; 
     _omega_1_1[19][0] = 0.84277; _omega_1_1[19][1] = 0.8446; _omega_1_1[19][2] = 0.8501; _omega_1_1[19][3] = 0.8592; _omega_1_1[19][4] = 0.8716; _omega_1_1[19][5] = 0.9053; _omega_1_1[19][6] = 0.9483; _omega_1_1[19][7] = 0.9984; 
     _omega_1_1[20][0] = 0.81287; _omega_1_1[20][1] = 0.8142; _omega_1_1[20][2] = 0.8183; _omega_1_1[20][3] = 0.8251; _omega_1_1[20][4] = 0.8344; _omega_1_1[20][5] = 0.8598; _omega_1_1[20][6] = 0.8927; _omega_1_1[20][7] = 0.9316; 
     _omega_1_1[21][0] = 0.78976; _omega_1_1[21][1] = 0.7908; _omega_1_1[21][2] = 0.7940; _omega_1_1[21][3] = 0.7993; _omega_1_1[21][4] = 0.8066; _omega_1_1[21][5] = 0.8265; _omega_1_1[21][6] = 0.8526; _omega_1_1[21][7] = 0.8836; 
     _omega_1_1[22][0] = 0.77111; _omega_1_1[22][1] = 0.7720; _omega_1_1[22][2] = 0.7745; _omega_1_1[22][3] = 0.7788; _omega_1_1[22][4] = 0.7846; _omega_1_1[22][5] = 0.8007; _omega_1_1[22][6] = 0.8219; _omega_1_1[22][7] = 0.8474; 
     _omega_1_1[23][0] = 0.75553; _omega_1_1[23][1] = 0.7562; _omega_1_1[23][2] = 0.7584; _omega_1_1[23][3] = 0.7619; _omega_1_1[23][4] = 0.7667; _omega_1_1[23][5] = 0.7800; _omega_1_1[23][6] = 0.7976; _omega_1_1[23][7] = 0.8189; 
     _omega_1_1[24][0] = 0.74220; _omega_1_1[24][1] = 0.7428; _omega_1_1[24][2] = 0.7446; _omega_1_1[24][3] = 0.7475; _omega_1_1[24][4] = 0.7515; _omega_1_1[24][5] = 0.7627; _omega_1_1[24][6] = 0.7776; _omega_1_1[24][7] = 0.7957; 
     _omega_1_1[25][0] = 0.72022; _omega_1_1[25][1] = 0.7206; _omega_1_1[25][2] = 0.7220; _omega_1_1[25][3] = 0.7241; _omega_1_1[25][4] = 0.7271; _omega_1_1[25][5] = 0.7354; _omega_1_1[25][6] = 0.7464; _omega_1_1[25][7] = 0.7600; 
     _omega_1_1[26][0] = 0.70254; _omega_1_1[26][1] = 0.7029; _omega_1_1[26][2] = 0.7039; _omega_1_1[26][3] = 0.7055; _omega_1_1[26][4] = 0.7078; _omega_1_1[26][5] = 0.7142; _omega_1_1[26][6] = 0.7228; _omega_1_1[26][7] = 0.7334; 
     _omega_1_1[27][0] = 0.68776; _omega_1_1[27][1] = 0.6880; _omega_1_1[27][2] = 0.6888; _omega_1_1[27][3] = 0.6901; _omega_1_1[27][4] = 0.6919; _omega_1_1[27][5] = 0.6970; _omega_1_1[27][6] = 0.7040; _omega_1_1[27][7] = 0.7125; 
     _omega_1_1[28][0] = 0.67510; _omega_1_1[28][1] = 0.6753; _omega_1_1[28][2] = 0.6760; _omega_1_1[28][3] = 0.6770; _omega_1_1[28][4] = 0.6785; _omega_1_1[28][5] = 0.6827; _omega_1_1[28][6] = 0.6884; _omega_1_1[28][7] = 0.6955; 
     _omega_1_1[29][0] = 0.66405; _omega_1_1[29][1] = 0.6642; _omega_1_1[29][2] = 0.6648; _omega_1_1[29][3] = 0.6657; _omega_1_1[29][4] = 0.6669; _omega_1_1[29][5] = 0.6704; _omega_1_1[29][6] = 0.6752; _omega_1_1[29][7] = 0.6811; 
     _omega_1_1[30][0] = 0.64136; _omega_1_1[30][1] = 0.6415; _omega_1_1[30][2] = 0.6418; _omega_1_1[30][3] = 0.6425; _omega_1_1[30][4] = 0.6433; _omega_1_1[30][5] = 0.6457; _omega_1_1[30][6] = 0.6490; _omega_1_1[30][7] = 0.6531; 
     _omega_1_1[31][0] = 0.62350; _omega_1_1[31][1] = 0.6236; _omega_1_1[31][2] = 0.6239; _omega_1_1[31][3] = 0.6243; _omega_1_1[31][4] = 0.6249; _omega_1_1[31][5] = 0.6267; _omega_1_1[31][6] = 0.6291; _omega_1_1[31][7] = 0.6321; 
     _omega_1_1[32][0] = 0.60882; _omega_1_1[32][1] = 0.6089; _omega_1_1[32][2] = 0.6091; _omega_1_1[32][3] = 0.6094; _omega_1_1[32][4] = 0.6099; _omega_1_1[32][5] = 0.6112; _omega_1_1[32][6] = 0.6131; _omega_1_1[32][7] = 0.6154; 
     _omega_1_1[33][0] = 0.59640; _omega_1_1[33][1] = 0.5964; _omega_1_1[33][2] = 0.5966; _omega_1_1[33][3] = 0.5969; _omega_1_1[33][4] = 0.5972; _omega_1_1[33][5] = 0.5983; _omega_1_1[33][6] = 0.5998; _omega_1_1[33][7] = 0.6017; 
     _omega_1_1[34][0] = 0.57626; _omega_1_1[34][1] = 0.5763; _omega_1_1[34][2] = 0.5764; _omega_1_1[34][3] = 0.5766; _omega_1_1[34][4] = 0.5768; _omega_1_1[34][5] = 0.5775; _omega_1_1[34][6] = 0.5785; _omega_1_1[34][7] = 0.5798; 
     _omega_1_1[35][0] = 0.54146; _omega_1_1[35][1] = 0.5415; _omega_1_1[35][2] = 0.5416; _omega_1_1[35][3] = 0.5416; _omega_1_1[35][4] = 0.5418; _omega_1_1[35][5] = 0.5421; _omega_1_1[35][6] = 0.5424; _omega_1_1[35][7] = 0.5429; 
     _omega_1_1[36][0] = 0.51803; _omega_1_1[36][1] = 0.5181; _omega_1_1[36][2] = 0.5182; _omega_1_1[36][3] = 0.5184; _omega_1_1[36][4] = 0.5184; _omega_1_1[36][5] = 0.5185; _omega_1_1[36][6] = 0.5186; _omega_1_1[36][7] = 0.5187; 

     //// <omega(2,2)*>
     _omega_2_2[0][0]  = 4.1005;  _omega_2_2[0][1]  = 4.266;  _omega_2_2[0][2]  = 4.833;  _omega_2_2[0][3]  = 5.742;  _omega_2_2[0][4]  = 6.729;  _omega_2_2[0][5]  = 8.624;  _omega_2_2[0][6]  = 10.34;  _omega_2_2[0][7]  = 11.89;
     _omega_2_2[1][0]  = 3.2626;  _omega_2_2[1][1]  = 3.305;  _omega_2_2[1][2]  = 3.516;  _omega_2_2[1][3]  = 3.914;  _omega_2_2[1][4]  = 4.433;  _omega_2_2[1][5]  = 5.570;  _omega_2_2[1][6]  = 6.637;  _omega_2_2[1][7]  = 7.618; 
     _omega_2_2[2][0]  = 2.8399;  _omega_2_2[2][1]  = 2.836;  _omega_2_2[2][2]  = 2.936;  _omega_2_2[2][3]  = 3.168;  _omega_2_2[2][4]  = 3.511;  _omega_2_2[2][5]  = 4.329;  _omega_2_2[2][6]  = 5.126;  _omega_2_2[2][7]  = 5.874; 
     _omega_2_2[3][0]  = 2.5310;  _omega_2_2[3][1]  = 2.522;  _omega_2_2[3][2]  = 2.586;  _omega_2_2[3][3]  = 2.749;  _omega_2_2[3][4]  = 3.004;  _omega_2_2[3][5]  = 3.640;  _omega_2_2[3][6]  = 4.282;  _omega_2_2[3][7]  = 4.895; 
     _omega_2_2[4][0]  = 2.2837;  _omega_2_2[4][1]  = 2.277;  _omega_2_2[4][2]  = 2.329;  _omega_2_2[4][3]  = 2.460;  _omega_2_2[4][4]  = 2.665;  _omega_2_2[4][5]  = 3.187;  _omega_2_2[4][6]  = 3.727;  _omega_2_2[4][7]  = 4.249; 
     _omega_2_2[5][0]  = 2.0838;  _omega_2_2[5][1]  = 2.081;  _omega_2_2[5][2]  = 2.130;  _omega_2_2[5][3]  = 2.243;  _omega_2_2[5][4]  = 2.417;  _omega_2_2[5][5]  = 2.862;  _omega_2_2[5][6]  = 3.329;  _omega_2_2[5][7]  = 3.786; 
     _omega_2_2[6][0]  = 1.9220;  _omega_2_2[6][1]  = 1.924;  _omega_2_2[6][2]  = 1.970;  _omega_2_2[6][3]  = 2.072;  _omega_2_2[6][4]  = 2.225;  _omega_2_2[6][5]  = 2.614;  _omega_2_2[6][6]  = 3.028;  _omega_2_2[6][7]  = 3.435; 
     _omega_2_2[7][0]  = 1.7902;  _omega_2_2[7][1]  = 1.795;  _omega_2_2[7][2]  = 1.840;  _omega_2_2[7][3]  = 1.934;  _omega_2_2[7][4]  = 2.070;  _omega_2_2[7][5]  = 2.417;  _omega_2_2[7][6]  = 2.788;  _omega_2_2[7][7]  = 3.156; 
     _omega_2_2[8][0]  = 1.6823;  _omega_2_2[8][1]  = 1.689;  _omega_2_2[8][2]  = 1.733;  _omega_2_2[8][3]  = 1.820;  _omega_2_2[8][4]  = 1.944;  _omega_2_2[8][5]  = 2.258;  _omega_2_2[8][6]  = 2.596;  _omega_2_2[8][7]  = 2.933; 
     _omega_2_2[9][0]  = 1.5929;  _omega_2_2[9][1]  = 1.601;  _omega_2_2[9][2]  = 1.644;  _omega_2_2[9][3]  = 1.725;  _omega_2_2[9][4]  = 1.838;  _omega_2_2[9][5]  = 2.124;  _omega_2_2[9][6]  = 2.435;  _omega_2_2[9][7]  = 2.746; 
     _omega_2_2[10][0] = 1.4551;  _omega_2_2[10][1] = 1.465;  _omega_2_2[10][2] = 1.504;  _omega_2_2[10][3] = 1.574;  _omega_2_2[10][4] = 1.670;  _omega_2_2[10][5] = 1.913;  _omega_2_2[10][6] = 2.181;  _omega_2_2[10][7] = 2.451; 
     _omega_2_2[11][0] = 1.3551;  _omega_2_2[11][1] = 1.365;  _omega_2_2[11][2] = 1.400;  _omega_2_2[11][3] = 1.461;  _omega_2_2[11][4] = 1.544;  _omega_2_2[11][5] = 1.754;  _omega_2_2[11][6] = 1.989;  _omega_2_2[11][7] = 2.228; 
     _omega_2_2[12][0] = 1.2800;  _omega_2_2[12][1] = 1.289;  _omega_2_2[12][2] = 1.321;  _omega_2_2[12][3] = 1.374;  _omega_2_2[12][4] = 1.447;  _omega_2_2[12][5] = 1.630;  _omega_2_2[12][6] = 1.838;  _omega_2_2[12][7] = 2.053; 
     _omega_2_2[13][0] = 1.2219;  _omega_2_2[13][1] = 1.231;  _omega_2_2[13][2] = 1.259;  _omega_2_2[13][3] = 1.306;  _omega_2_2[13][4] = 1.370;  _omega_2_2[13][5] = 1.532;  _omega_2_2[13][6] = 1.718;  _omega_2_2[13][7] = 1.912; 
     _omega_2_2[14][0] = 1.1757;  _omega_2_2[14][1] = 1.184;  _omega_2_2[14][2] = 1.209;  _omega_2_2[14][3] = 1.251;  _omega_2_2[14][4] = 1.307;  _omega_2_2[14][5] = 1.451;  _omega_2_2[14][6] = 1.618;  _omega_2_2[14][7] = 1.795; 
     _omega_2_2[15][0] = 1.0933;  _omega_2_2[15][1] = 1.100;  _omega_2_2[15][2] = 1.119;  _omega_2_2[15][3] = 1.150;  _omega_2_2[15][4] = 1.193;  _omega_2_2[15][5] = 1.304;  _omega_2_2[15][6] = 1.435;  _omega_2_2[15][7] = 1.578; 
     _omega_2_2[16][0] = 1.0388;  _omega_2_2[16][1] = 1.044;  _omega_2_2[16][2] = 1.059;  _omega_2_2[16][3] = 1.083;  _omega_2_2[16][4] = 1.117;  _omega_2_2[16][5] = 1.204;  _omega_2_2[16][6] = 1.310;  _omega_2_2[16][7] = 1.428; 
     _omega_2_2[17][0] = 0.99938; _omega_2_2[17][1] = 1.004;  _omega_2_2[17][2] = 1.016;  _omega_2_2[17][3] = 1.035;  _omega_2_2[17][4] = 1.062;  _omega_2_2[17][5] = 1.133;  _omega_2_2[17][6] = 1.220;  _omega_2_2[17][7] = 1.319; 
     _omega_2_2[18][0] = 0.96988; _omega_2_2[18][1] = 0.9732; _omega_2_2[18][2] = 0.9830; _omega_2_2[18][3] = 0.9991; _omega_2_2[18][4] = 1.021;  _omega_2_2[18][5] = 1.079;  _omega_2_2[18][6] = 1.153;  _omega_2_2[18][7] = 1.236; 
     _omega_2_2[19][0] = 0.92676; _omega_2_2[19][1] = 0.9291; _omega_2_2[19][2] = 0.9360; _omega_2_2[19][3] = 0.9473; _omega_2_2[19][4] = 0.9628; _omega_2_2[19][5] = 1.005;  _omega_2_2[19][6] = 1.058;  _omega_2_2[19][7] = 1.121; 
     _omega_2_2[20][0] = 0.89616; _omega_2_2[20][1] = 0.8979; _omega_2_2[20][2] = 0.9030; _omega_2_2[20][3] = 0.9114; _omega_2_2[20][4] = 0.9230; _omega_2_2[20][5] = 0.9545; _omega_2_2[20][6] = 0.9955; _omega_2_2[20][7] = 1.044; 
     _omega_2_2[21][0] = 0.87272; _omega_2_2[21][1] = 0.8741; _omega_2_2[21][2] = 0.8780; _omega_2_2[21][3] = 0.8845; _omega_2_2[21][4] = 0.8935; _omega_2_2[21][5] = 0.9181; _omega_2_2[21][6] = 0.9505; _omega_2_2[21][7] = 0.9893; 
     _omega_2_2[22][0] = 0.85379; _omega_2_2[22][1] = 0.8549; _omega_2_2[22][2] = 0.8580; _omega_2_2[22][3] = 0.8632; _omega_2_2[22][4] = 0.8703; _omega_2_2[22][5] = 0.8901; _omega_2_2[22][6] = 0.9164; _omega_2_2[22][7] = 0.9482; 
     _omega_2_2[23][0] = 0.83795; _omega_2_2[23][1] = 0.8388; _omega_2_2[23][2] = 0.8414; _omega_2_2[23][3] = 0.8456; _omega_2_2[23][4] = 0.8515; _omega_2_2[23][5] = 0.8678; _omega_2_2[23][6] = 0.8895; _omega_2_2[23][7] = 0.9160; 
     _omega_2_2[24][0] = 0.82435; _omega_2_2[24][1] = 0.8251; _omega_2_2[24][2] = 0.8273; _omega_2_2[24][3] = 0.8308; _omega_2_2[24][4] = 0.8356; _omega_2_2[24][5] = 0.8493; _omega_2_2[24][6] = 0.8676; _omega_2_2[24][7] = 0.8901; 
     _omega_2_2[25][0] = 0.80184; _omega_2_2[25][1] = 0.8024; _omega_2_2[25][2] = 0.8039; _omega_2_2[25][3] = 0.8065; _omega_2_2[25][4] = 0.8101; _omega_2_2[25][5] = 0.8201; _omega_2_2[25][6] = 0.8337; _omega_2_2[25][7] = 0.8506; 
     _omega_2_2[26][0] = 0.78363; _omega_2_2[26][1] = 0.7840; _omega_2_2[26][2] = 0.7852; _omega_2_2[26][3] = 0.7872; _omega_2_2[26][4] = 0.7899; _omega_2_2[26][5] = 0.7976; _omega_2_2[26][6] = 0.8081; _omega_2_2[26][7] = 0.8212; 
     _omega_2_2[27][0] = 0.76834; _omega_2_2[27][1] = 0.7687; _omega_2_2[27][2] = 0.7696; _omega_2_2[27][3] = 0.7712; _omega_2_2[27][4] = 0.7733; _omega_2_2[27][5] = 0.7794; _omega_2_2[27][6] = 0.7878; _omega_2_2[27][7] = 0.7983; 
     _omega_2_2[28][0] = 0.75518; _omega_2_2[28][1] = 0.7554; _omega_2_2[28][2] = 0.7562; _omega_2_2[28][3] = 0.7575; _omega_2_2[28][4] = 0.7592; _omega_2_2[28][5] = 0.7642; _omega_2_2[28][6] = 0.7711; _omega_2_2[28][7] = 0.7797; 
     _omega_2_2[29][0] = 0.74364; _omega_2_2[29][1] = 0.7438; _omega_2_2[29][2] = 0.7445; _omega_2_2[29][3] = 0.7455; _omega_2_2[29][4] = 0.7470; _omega_2_2[29][5] = 0.7512; _omega_2_2[29][6] = 0.7569; _omega_2_2[29][7] = 0.7642; 
     _omega_2_2[30][0] = 0.71982; _omega_2_2[30][1] = 0.7200; _omega_2_2[30][2] = 0.7204; _omega_2_2[30][3] = 0.7211; _omega_2_2[30][4] = 0.7221; _omega_2_2[30][5] = 0.7250; _omega_2_2[30][6] = 0.7289; _omega_2_2[30][7] = 0.7339; 
     _omega_2_2[31][0] = 0.70097; _omega_2_2[31][1] = 0.7011; _omega_2_2[31][2] = 0.7014; _omega_2_2[31][3] = 0.7019; _omega_2_2[31][4] = 0.7026; _omega_2_2[31][5] = 0.7047; _omega_2_2[31][6] = 0.7076; _omega_2_2[31][7] = 0.7112; 
     _omega_2_2[32][0] = 0.68545; _omega_2_2[32][1] = 0.6855; _omega_2_2[32][2] = 0.6858; _omega_2_2[32][3] = 0.6861; _omega_2_2[32][4] = 0.6867; _omega_2_2[32][5] = 0.6883; _omega_2_2[32][6] = 0.6905; _omega_2_2[32][7] = 0.6932; 
     _omega_2_2[33][0] = 0.67232; _omega_2_2[33][1] = 0.6724; _omega_2_2[33][2] = 0.6726; _omega_2_2[33][3] = 0.6728; _omega_2_2[33][4] = 0.6733; _omega_2_2[33][5] = 0.6745; _omega_2_2[33][6] = 0.6762; _omega_2_2[33][7] = 0.6784; 
     _omega_2_2[34][0] = 0.65099; _omega_2_2[34][1] = 0.6510; _omega_2_2[34][2] = 0.6512; _omega_2_2[34][3] = 0.6513; _omega_2_2[34][4] = 0.6516; _omega_2_2[34][5] = 0.6524; _omega_2_2[34][6] = 0.6534; _omega_2_2[34][7] = 0.6546; 
     _omega_2_2[35][0] = 0.61397; _omega_2_2[35][1] = 0.6141; _omega_2_2[35][2] = 0.6143; _omega_2_2[35][3] = 0.6145; _omega_2_2[35][4] = 0.6147; _omega_2_2[35][5] = 0.6148; _omega_2_2[35][6] = 0.6148; _omega_2_2[35][7] = 0.6147; 
     _omega_2_2[36][0] = 0.58870; _omega_2_2[36][1] = 0.5889; _omega_2_2[36][2] = 0.5894; _omega_2_2[36][3] = 0.5900; _omega_2_2[36][4] = 0.5903; _omega_2_2[36][5] = 0.5901; _omega_2_2[36][6] = 0.5895; _omega_2_2[36][7] = 0.5885; 
  }


} //end namespace Antioch

#endif
