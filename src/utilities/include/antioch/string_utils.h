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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef ANTIOCH_STRING_UTILS_H
#define ANTIOCH_STRING_UTILS_H

// Antioch
#include "antioch/antioch_asserts.h"
// some enum
#include "antioch/kinetics_enum.h"
#include "antioch/reaction_enum.h"

// C++
#include <iostream>
#include <string>
#include <vector>
#include <sstream>

namespace Antioch
{
  //! All characters in delimiter will be treated as a delimiter
  void split_string( const std::string& input,
                     const std::string& delimiter,
                     std::vector<std::string>& results );


  template <typename T>
  inline
  T string_to_T(const std::string& input)
  {
    std::istringstream converter(input);
    T returnval;
    converter >> returnval;
    if (converter.fail())
      antioch_error();
    return returnval;
  }

  template <typename Type>
  inline
  std::pair<std::string, Type> split_string_on_colon(const std::string &token)
  {
    std::pair<std::string, Type> ret = std::make_pair(std::string(), 0);
    std::string::size_type colon_position = token.find(":");
    antioch_assert (colon_position != std::string::npos);
    ret.first  = token.substr(0, colon_position);
    ret.second = string_to_T<Type>(token.substr(colon_position + 1));
    return ret;
  }

  /*!
    Taken from FIN-S for XML parsing.
   */
  inline
  int SplitString(const std::string& input,
		  const std::string& delimiter,
		  std::vector<std::string>& results,
		  bool includeEmpties = true)
  {
    using std::vector;
    using std::string;

    int iPos = 0;
    int newPos = -1;
    int sizeS2 = (int)delimiter.size();
    int isize = (int)input.size();

    if(
       ( isize == 0 )
       ||
       ( sizeS2 == 0 )
	)
      {
	return 0;
      }

    vector<int> positions;

    newPos = input.find (delimiter, 0);

    if( newPos < 0 )
      {
	return 0;
      }

    int numFound = 0;

    while( newPos >= iPos )
      {
	numFound++;
	positions.push_back(newPos);
	iPos = newPos;
	newPos = input.find (delimiter, iPos+sizeS2);
      }

    if( numFound == 0 )
      {
	return 0;
      }

    for( int i=0; i <= static_cast<int>(positions.size()); ++i )
      {
	string s("");
	if( i == 0 )
	  {
	    s = input.substr( i, positions[i] );
	  }
	else
	  {
	    int offset = positions[i-1] + sizeS2;
	    if( offset < isize )
	      {
		if( i == static_cast<int>(positions.size()) )
		  {
		    s = input.substr(offset);
		  }
		else if( i > 0 )
		  {
		    s = input.substr( positions[i-1] + sizeS2,
				      positions[i] - positions[i-1] - sizeS2 );
		  }
	      }
	  }
	if( includeEmpties || ( s.size() > 0 ) )
	  {
	    results.push_back(s);
	  }
      }
    return numFound;
  }

  /*!
      adapted getline, never believe ascii file for the
      formatting of end-of-line.
      end-of-line triggered by \n or \r
        - Windows     \n\r
        - Unix/Linux  \n
        - Mac         \r
   */
  inline
  int ascii_getline(std::istream & buf, std::string & line)
  {
     char c('a');
     line.clear();
     while(!buf.eof())
     {
        c = buf.get();
        if(c == '\n' || c == '\r')break;
        line += c;
     }
     char n = buf.peek();

     /* never trust ascii files, they may come from
        Windows, suppodedly \n\r, but let's not
        underestimate Windows's viciousness
      */
     if((c == '\n' && n == '\r') ||
        (n == '\n' && c == '\r'))c = buf.get();

     return buf.good();
  }

  inline
  KineticsModel::Parameters string_to_kin_enum(const std::string & str)
  {
// kinetics
      if(str == "A")
      {
        return KineticsModel::Parameters::A;
      }else if(str == "E")
      {
        return KineticsModel::Parameters::E;
      }else if(str == "B")
      {
        return KineticsModel::Parameters::B;
      }else if(str == "D")
      {
        return KineticsModel::Parameters::D;
      }else if(str == "Tref")
      {
        return KineticsModel::Parameters::T_REF;
      }else if(str == "Rscale")
      {
        return KineticsModel::Parameters::R_SCALE;
      }else if(str == "sigma")
      {
        return KineticsModel::Parameters::SIGMA;
      }else if(str == "lambda")
      {
        return KineticsModel::Parameters::LAMBDA;
      }else if(str == "0")
      {
        return KineticsModel::Parameters::LOW_PRESSURE;
      }else if(str == "inf")
      {
        return KineticsModel::Parameters::HIGH_PRESSURE;
      }else
      {
        return KineticsModel::Parameters::NOT_FOUND;
      }
   }

  inline
  ReactionType::Parameters string_to_chem_enum(const std::string & str)
  {
// chemical
      if(str == "efficiencies")
      {
        return ReactionType::Parameters::EFFICIENCIES;
      }else if(str == "alpha")
      {
        return ReactionType::Parameters::TROE_ALPHA;
      }else if(str == "T1")
      {
        return ReactionType::Parameters::TROE_T1;
      }else if(str == "T2")
      {
        return ReactionType::Parameters::TROE_T2;
      }else if(str == "T3")
      {
        return ReactionType::Parameters::TROE_T3;
      }else
      {
        return ReactionType::Parameters::NOT_FOUND;
      }
  }

} // end namespace Antioch

#endif // ANTIOCH_STRING_UTILS_H
