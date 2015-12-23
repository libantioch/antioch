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

//Antioch

//C++
#include <string>
#include <map>

namespace Antioch
{

/*ChemKin spec*/

  class ChemKinDefinitions
  {
    public:
    ChemKinDefinitions();
    ~ChemKinDefinitions(){}

    enum Delim{ ERROR = -1,
                PLUS,
                IRREVERSIBLE,
                REVERSIBLE,
                REVERSIBLE_ALT
              };

    enum Symbol{
                 TB = 0,
                 FAL,
                 PHOTO,
                 ELECTRO
             };

    const std::map<Delim,std::string> & delim()   const;

    const std::map<Symbol,std::string> & symbol() const;

    const std::string & reversible()              const;

    const std::string & duplicate()               const;

    const std::string & end_tag()                 const;

    const std::string & comment()                 const;

    const std::string & parser()                  const;

    bool is_comment(const char & c)               const;

    bool is_equation_delimiter(const std::string & test) const;

    Delim equation_delimiter(const std::string & test) const;

    private:
      const std::string            _reversible;
      const std::string            _duplicate;
      const std::string            _end_tag;
      const std::string            _comment;
      const std::string            _parser;
      std::map<Delim,std::string>  _delim;
      std::map<Symbol,std::string> _symbol;

  };

  inline
  ChemKinDefinitions::ChemKinDefinitions():
     _reversible("REV"),
     _duplicate("DUP"),
     _end_tag("END"),
     _comment("!"),
     _parser("/")
  {
    _delim[PLUS]           = "+";
    _delim[REVERSIBLE]     = "=";
    _delim[REVERSIBLE_ALT] = "<=>";
    _delim[IRREVERSIBLE]   = "=>";

    _symbol[TB]      = "M";
    _symbol[FAL]     = "(+M)";
    _symbol[PHOTO]   = "HV";
    _symbol[ELECTRO] = "E";
    _symbol[FORD]    = "FORD";
    _symbol[RORD]    = "RORD";
  }

  inline
  const std::map<ChemKinDefinitions::Delim,std::string> & ChemKinDefinitions::delim()   const 
  {
    return _delim;
  }

  inline
  const std::map<ChemKinDefinitions::Symbol,std::string> & ChemKinDefinitions::symbol() const 
  {
    return _symbol;
  }

  inline
  const std::string & ChemKinDefinitions::reversible() const 
  {
    return _reversible;
  }

  inline
  const std::string & ChemKinDefinitions::duplicate() const 
  {
    return _duplicate;
  }

  inline
  const std::string & ChemKinDefinitions::end_tag() const 
  {
    return _end_tag;
  }

  inline
  const std::string & ChemKinDefinitions::comment() const 
  {
    return _comment;
  }

  inline
  const std::string & ChemKinDefinitions::parser() const 
  {
    return _parser;
  }

  inline
  bool ChemKinDefinitions::is_comment(const char & c) const 
  {
    return (c == _comment[0]);
  }

  inline
  bool ChemKinDefinitions::is_equation_delimiter(const std::string & test) const 
  {
    return (test == _delim.at(REVERSIBLE)     || 
            test == _delim.at(REVERSIBLE_ALT) || 
            test == _delim.at(IRREVERSIBLE)
           );
  }

  inline
  ChemKinDefinitions::Delim ChemKinDefinitions::equation_delimiter(const std::string & test) const 
  {
    if(test.find(_delim.at(REVERSIBLE_ALT)) != std::string::npos)
    {
       return REVERSIBLE_ALT;
    }else if(test.find(_delim.at(IRREVERSIBLE)) != std::string::npos)
    {
      return IRREVERSIBLE;
    }else if(test.find(_delim.at(REVERSIBLE)) != std::string::npos)
    {
      return REVERSIBLE;
    }else
    {
      return ERROR;
    }
  }

}
