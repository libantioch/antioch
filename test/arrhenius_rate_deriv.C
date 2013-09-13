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

template <typename Scalar>
int tester()
{
  using std::abs;
  using std::exp;
  using std::pow;
  using MetaPhysicL::DualNumber;

  const Scalar Cf = 1.4;
  const Scalar eta = 1.2;
  const Scalar Ea = 5.0;

  Antioch::ArrheniusRate<Scalar> arrhenius_rate(Cf,eta,Ea);
  Antioch::ArrheniusRate<DualNumber<Scalar> > arrhenius_T_sensitivity(Cf,eta,Ea);

  const Scalar T = 1500.1;
  const DualNumber<Scalar> T_DN(T, 1);
  
  int return_flag = 0;

  Scalar rate, deriveRate;

  arrhenius_rate.rate_and_derivative(T,rate,deriveRate);
  DualNumber<Scalar> rate_DN_T = arrhenius_T_sensitivity(T_DN);

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;

  if( abs( (deriveRate - rate_DN_T.derivatives())/deriveRate) > tol )
    {
      std::cout << "Error: Mismatch in rate derivatives." << std::endl
                << std::setprecision(25)
		<< "deriveRate(T) = " << deriveRate << std::endl
		<< "rate_DN_T' =    " << rate_DN_T.derivatives() << std::endl;

      return_flag = 1;
    }

  return return_flag;
}
#endif

int main()
{
#ifdef ANTIOCH_HAVE_METAPHYSICL
  return (tester<double>() ||
          tester<long double>() ||
          tester<float>());
#else
  return 77; // automake code for a skipped test
#endif
}
