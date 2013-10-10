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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "antioch_config.h"

#ifdef ANTIOCH_HAVE_METAPHYSICL
#include "metaphysicl/dualnumber.h"

#include "antioch/arrhenius_rate.h"

#include <iomanip>
#include <string>

template <typename Scalar>
int tester(const std::string& testname)
{
  using std::abs;
  using std::exp;
  using std::pow;
  using MetaPhysicL::DualNumber;

  int return_flag = 0;

  const Scalar Cf = 1.4;
  const Scalar eta = 1.2;
  const Scalar Ea = 5.0;

  const DualNumber<Scalar> Ea_DN(Ea, 1);

  Antioch::ArrheniusRate<Scalar> arrhenius_rate(Cf,eta,Ea);

  const Scalar T = 1500.1;

  Scalar rate, deriveRate;

  arrhenius_rate.rate_and_derivative(T,rate,deriveRate);


  Antioch::ArrheniusRate<DualNumber<Scalar> >
    arrhenius_T_sensitivity(Cf,eta,Ea);

  const DualNumber<Scalar> T_DN(T, 1);

  DualNumber<Scalar> rate_DN_T = arrhenius_T_sensitivity(T_DN);

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;

  if( abs( (deriveRate - rate_DN_T.derivatives())/deriveRate) > tol )
    {
      std::cout << "Error: Mismatch in " << testname << " rate derivatives." << std::endl
                << std::setprecision(25)
		<< "deriveRate(T) = " << deriveRate << std::endl
		<< "rate_DN_T' =    " << rate_DN_T.derivatives() << std::endl;

      return_flag = 1;
    }


  Antioch::ArrheniusRate<DualNumber<Scalar> >
    arrhenius_Ea_sensitivity(Cf,eta,Ea_DN);

  const Scalar perturb =
          std::pow(std::numeric_limits<Scalar>::epsilon(), (1./3.)) * 100;

  Antioch::ArrheniusRate<Scalar> arrhenius_rate_Eaplus(Cf,eta,Ea + perturb);
  Antioch::ArrheniusRate<Scalar> arrhenius_rate_Eaminus(Cf,eta,Ea - perturb);

  // FIXME - we need to upgrade the return type of templated functions
  // in the case where the input StateType is *weaker* than CoeffType
  DualNumber<Scalar> rate_DN_Ea =
          arrhenius_Ea_sensitivity(DualNumber<Scalar>(T) );

  Scalar dR_dEa = (arrhenius_rate_Eaplus(T) -
                   arrhenius_rate_Eaminus(T)) / (2 * perturb);

  if( abs( (dR_dEa - rate_DN_Ea.derivatives())/dR_dEa) > 
      perturb * perturb * 100)
    {
      std::cout << "Error: Mismatch in " << testname << " rate derivatives." << std::endl
                << std::setprecision(25)
		<< "perturb =      " << perturb << std::endl
		<< "dR_dEa =      " << dR_dEa << std::endl
		<< "rate_DN_Ea' = " << rate_DN_Ea.derivatives() << std::endl;

      return_flag = 1;
    }

  return return_flag;
}
#endif

int main()
{
#ifdef ANTIOCH_HAVE_METAPHYSICL
  return (tester<double>("double") ||
          tester<long double>("long double") ||
          tester<float>("float"));
#else
  return 77; // automake code for a skipped test
#endif
}
