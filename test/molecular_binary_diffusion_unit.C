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

// C++
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

// Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/transport_species.h"
#include "antioch/molecular_binary_diffusion.h"
#include "antioch/vector_utils.h"

template <typename Scalar>
int test_diff( const Scalar dij, const Scalar dij_exact, const Scalar tol, const std::string & words )
{
  using std::abs;

  int return_flag = 0;

  const double rel_error = abs( (dij - dij_exact)/dij_exact);

  if( rel_error  > tol )
    {
      std::cerr << std::setprecision(15) << std::scientific;
      std::cerr << "Error: Mismatch in bimolecular coefficient of " << words << std::endl
		<< "Dij(T)    = " << dij << std::endl
		<< "Dij_exact = " << dij_exact << std::endl
		<< "rel_error = " << rel_error << std::endl
		<< "tol = " << tol << std::endl;
      return_flag = 1;
    }

  return return_flag;
}

template <typename Scalar>
int tester()
{
/* from default data
# Species, eps/kB (K), sigma (ang), alpha (D), alpha (ang^3), Zrot@298 K
N2                     97.530     3.621     0.000     1.760     4.000
CH4                   141.400     3.746     0.000     2.600    13.000
H2O                   572.400     2.605     1.844     0.000     4.000
*/

//  N2
   const Scalar N2_LJ_eps(97.530L);
   const Scalar N2_LJ_depth(3.621L);
   const Scalar N2_dipole(0.L);
   const Scalar N2_polar(1.760L);
   const Scalar N2_Zrot(4.L);
   const Scalar N2_mass(28.016e-3L);
// CH4
   const Scalar CH4_LJ_eps(141.400L);
   const Scalar CH4_LJ_depth(3.746L);
   const Scalar CH4_dipole(0.L);
   const Scalar CH4_polar(2.6L);
   const Scalar CH4_Zrot(13.L);
   const Scalar CH4_mass(16.043e-3L);
// H2O
   const Scalar H2O_LJ_eps(572.4L);
   const Scalar H2O_LJ_depth(2.605L);
   const Scalar H2O_dipole(1.844L);
   const Scalar H2O_polar(0.L);
   const Scalar H2O_Zrot(4.L);
   const Scalar H2O_mass(18.016e-3L);
  Antioch::TransportSpecies<Scalar> N2(0,N2_LJ_eps,N2_LJ_depth,N2_dipole,N2_polar,N2_Zrot,N2_mass),
                                    CH4(1,CH4_LJ_eps,CH4_LJ_depth,CH4_dipole,CH4_polar,CH4_Zrot,CH4_mass),
                                    H2O(2,H2O_LJ_eps,H2O_LJ_depth,H2O_dipole,H2O_polar,H2O_Zrot,H2O_mass);


  Antioch::MolecularBinaryDiffusion<Scalar,Antioch::GSLSpliner> 
                D00(N2,N2),  D01(N2,CH4),  D02(N2,H2O),
                D10(CH4,N2), D11(CH4,CH4), D12(CH4,H2O),
                D20(H2O,N2), D21(H2O,CH4), D22(H2O,H2O);

  const Scalar T    = 1500.1;
  const Scalar cTot = 5e-7;

  // bc gives
  const Scalar D00_exact_times_interp = 3.575629059282712203350417538824313603e3;
  const Scalar D11_exact_times_interp = 4.415035849326582722446637772820322872e3;
  const Scalar D22_exact_times_interp = 8.615250909281137767894964155009068175e3;
  const Scalar D10_exact_times_interp = 4.048999169845892614961423528844221310e3;
  const Scalar D20_exact_times_interp = 5.468107000169991297723144486054211050e3;
  const Scalar D12_exact_times_interp = 5.973341001459783642751059656941311432e3;

  int return_flag = 0;

  const Scalar tol = (std::numeric_limits<Scalar>::epsilon() * 10. < 2e-16)?2e-16:
                                                                            std::numeric_limits<Scalar>::epsilon() * 10.;

// symetric consistency
  std::cout << "***  Testing symetry ...";
            
  return_flag = test_diff( D10(T,cTot) , D01(T,cTot), tol, "N2 - CH4 symetry" )  || return_flag;
  return_flag = test_diff( D12(T,cTot) , D21(T,cTot), tol, "CH4 - H2O symetry" ) || return_flag;
  return_flag = test_diff( D20(T,cTot) , D02(T,cTot), tol, "H2O - N2 symetry" )  || return_flag;

  (return_flag)?std::cout << " ...failed\n":std::cout << " ...passed\n";
  std::cout << "*****************************\n" << std::endl;

// values
  return_flag = test_diff( D00(T,cTot) * D00.Stockmayer(T), D00_exact_times_interp, tol, "N2 - N2" )   || return_flag;
  return_flag = test_diff( D11(T,cTot) * D11.Stockmayer(T), D11_exact_times_interp, tol, "CH4 - CH4" ) || return_flag;
  return_flag = test_diff( D22(T,cTot) * D22.Stockmayer(T), D22_exact_times_interp, tol, "H2O - H2O" ) || return_flag;
  return_flag = test_diff( D10(T,cTot) * D10.Stockmayer(T), D10_exact_times_interp, tol, "N2 - CH4" )  || return_flag;
  return_flag = test_diff( D12(T,cTot) * D12.Stockmayer(T), D12_exact_times_interp, tol, "CH4 - H2O" ) || return_flag;
  return_flag = test_diff( D20(T,cTot) * D20.Stockmayer(T), D20_exact_times_interp, tol, "N2 - H2O" )  || return_flag;

  return return_flag;
}

int main()
{
  return tester<double>()    ||
         tester<long double>() ||
         tester<float>();
};
