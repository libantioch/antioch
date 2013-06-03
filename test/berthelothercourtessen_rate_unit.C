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
// $Id: arrhenius_rate_unit.C 38747 2013-04-17 23:26:39Z splessis $
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "antioch/berthelothercourtessen_rate.h"

template <typename Scalar>
int tester()
{
  using std::abs;
  using std::exp;
  using std::pow;

  const Scalar Cf  = 1.4;
  const Scalar eta = 1.2;
  const Scalar D   = -5.0;

  Antioch::BerthelotHercourtEssenRate<Scalar> berthelothercourtessen_rate(Cf,eta,D);

  const Scalar T = 1500.1;
  
  const Scalar rate_exact = Cf*pow(T,eta)*exp(D*T);

  int return_flag = 0;

  Scalar rate = berthelothercourtessen_rate(T);

  const Scalar tol = 1.0e-15;

  if( abs( (rate - rate_exact)/rate_exact ) > tol )
    {
      std::cout << "Error: Mismatch in rate values." << std::endl
		<< "rate(T) = " << rate << std::endl
		<< "rate_exact = " << rate_exact << std::endl;

      return_flag = 1;
    }

  std::cout << "Berthelot-Hercourt-Essen rate: " << berthelothercourtessen_rate << std::endl;

  return return_flag;
}

int main()
{
  return (tester<double>() ||
          tester<long double>() ||
          tester<float>());
}
