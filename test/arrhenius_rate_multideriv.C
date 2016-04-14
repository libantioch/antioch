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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "antioch_config.h"

#ifdef ANTIOCH_HAVE_METAPHYSICL
#include "metaphysicl/dualnumber.h"
#include "metaphysicl/dynamicsparsenumberarray.h"

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

  typedef MetaPhysicL::DynamicSparseNumberArray<Scalar, int> DNSA;
  typedef DualNumber<Scalar, DNSA> DN;

  int return_flag = 0;

  const Scalar Cf = 1.4;
  const Scalar Ea = 5.0;
  const Scalar rscale = Antioch::Constants::R_universal<Scalar>();

  // Unit vectors for independent variable derivatives
  DNSA unit_Cf, unit_Ea, unit_rscale;
  unit_Cf.resize(1);
  unit_Cf.raw_index(0) = 0;
  unit_Cf.raw_at(0) = 1;
  unit_Ea.resize(1);
  unit_Ea.raw_index(0) = 1;
  unit_Ea.raw_at(0) = 1;
  unit_rscale.resize(1);
  unit_rscale.raw_index(0) = 2;
  unit_rscale.raw_at(0) = 1;

  const DN Cf_DN(Cf, unit_Cf);
  const DN Ea_DN(Ea, unit_Ea);
  const DN rscale_DN(rscale, unit_rscale);

  const DN partsum_DN = Ea_DN + Cf_DN;
  const DN fullsum_DN = Ea_DN + Cf_DN + rscale_DN;

  const Scalar T = 1500.1;

  const Antioch::ArrheniusRate<Scalar> arrhenius_rate(Cf, Ea, rscale);

  const Scalar rate = arrhenius_rate(T);

  const Antioch::ArrheniusRate<DN>
    arrhenius_param_sensitivity(Cf_DN, Ea_DN, rscale_DN);

  const Scalar perturb =
          std::pow(std::numeric_limits<Scalar>::epsilon(), (1./3.)) * 100;

  const Antioch::ArrheniusRate<Scalar>
    arrhenius_rate_Cf_plus(Cf + perturb, Ea, rscale);
  const Antioch::ArrheniusRate<Scalar>
    arrhenius_rate_Cf_minus(Cf - perturb, Ea, rscale);

  const Antioch::ArrheniusRate<Scalar>
    arrhenius_rate_Ea_plus(Cf, Ea + perturb, rscale);
  const Antioch::ArrheniusRate<Scalar>
    arrhenius_rate_Ea_minus(Cf, Ea - perturb, rscale);

  const Antioch::ArrheniusRate<Scalar>
    arrhenius_rate_rscale_plus(Cf, Ea, rscale + perturb);
  const Antioch::ArrheniusRate<Scalar>
    arrhenius_rate_rscale_minus(Cf, Ea, rscale - perturb);

  // FIXME - we need to upgrade the return type of templated functions
  // in the case where the input StateType is *weaker* than CoeffType
  const DN rate_DN = arrhenius_param_sensitivity(DN(T));

  const Scalar dR_dCf = (arrhenius_rate_Cf_plus(T) -
                         arrhenius_rate_Cf_minus(T)) / (2 * perturb);

  const Scalar dR_dEa = (arrhenius_rate_Ea_plus(T) -
                         arrhenius_rate_Ea_minus(T)) / (2 * perturb);

  const Scalar dR_drscale = (arrhenius_rate_rscale_plus(T) -
                             arrhenius_rate_rscale_minus(T)) / (2 * perturb);

  if( abs( (rate - rate_DN.value())/rate) > 
      std::numeric_limits<Scalar>::epsilon() * 100)
    {
      std::cout << "Error: Mismatch in " << testname << " rate." << std::endl
                << std::setprecision(25)
		<< "rate =    " << rate << std::endl
		<< "rate_DN = " << rate_DN.value() << std::endl;

      return_flag = 1;
    }


  if( abs( (dR_dCf - rate_DN.derivatives()[0])/dR_dCf) > 
      perturb * perturb * 100)
    {
      std::cout << "Error: Mismatch in " << testname << " rate derivatives." << std::endl
                << std::setprecision(25)
		<< "perturb =      " << perturb << std::endl
		<< "dR_dCf =      " << dR_dCf << std::endl
		<< "drate_DN/dCf' = " << rate_DN.derivatives()[0] << std::endl;

      return_flag = 1;
    }


  if( abs( (dR_dEa - rate_DN.derivatives()[1])/dR_dEa) > 
      perturb * perturb * 100)
    {
      std::cout << "Error: Mismatch in " << testname << " rate derivatives." << std::endl
                << std::setprecision(25)
		<< "perturb =      " << perturb << std::endl
		<< "dR_dEa =      " << dR_dEa << std::endl
		<< "drate_DN/dEa' = " << rate_DN.derivatives()[1] << std::endl;

      return_flag = 1;
    }


  if( abs( (dR_drscale - rate_DN.derivatives()[2])/dR_drscale) > 
      perturb * perturb * 100)
    {
      std::cout << "Error: Mismatch in " << testname << " rate derivatives." << std::endl
                << std::setprecision(25)
		<< "perturb =      " << perturb << std::endl
		<< "dR_drscale =      " << dR_drscale << std::endl
		<< "drate_DN/drscale' = " << rate_DN.derivatives()[1] << std::endl;

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
