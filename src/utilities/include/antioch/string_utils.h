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

// C++
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib> // atoi

namespace Antioch
{
  /*!
    Split on colon, and return name, int value pair.
    Taken from FIN-S for XML parsing.
   */
  inline
  std::pair<std::string, int> split_string_int_on_colon(const std::string &token)
  {
    std::pair<std::string, int> ret = std::make_pair(std::string(), 0);
    std::string::size_type colon_position = token.find(":");
    antioch_assert (colon_position != std::string::npos);
    ret.first  = token.substr(0, colon_position);
    ret.second = std::atoi(token.substr(colon_position + 1).c_str());
    return ret;
  }


  /*!
    Split on colon, and return name, double value pair.
    Taken from FIN-S for XML parsing.
   */
  inline
  std::pair<std::string, double> split_string_double_on_colon(const std::string &token)
  {
    std::pair<std::string, double> ret = std::make_pair(std::string(), 0.0);
    std::string::size_type colon_position = token.find(":");
    antioch_assert (colon_position != std::string::npos);
    ret.first  = token.substr(0, colon_position);
    ret.second = std::atof(token.substr(colon_position + 1).c_str());
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


} // end namespace Antioch

#endif // ANTIOCH_STRING_UTILS_H
