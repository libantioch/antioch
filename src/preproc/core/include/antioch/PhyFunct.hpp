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
#ifndef _PHYSICAL_FUNCTIONS_
#define _PHYSICAL_FUNCTIONS_

#include "ParameterPhy.hpp"
#include <string>

namespace Antioch
{

/*
 * \class PhyFunct
 * \brief Regroups two physical parameters and a name.
 * 
 * This is useful for the class PhyHarmFunct. It is
 * useful for a user that may want to classify two
 * parameters as being one entity, typically to store
 * a function \f$ y = f(x)\f$, one parameter being \f$y\f$,
 * the other \f$x\f$.
 *
 * It is important to notice that this is but a structure, there
 * is absolutely no testing of anything, that the number of
 * stored values for \f$x\f$ and \f$y\f$ corresponds for instance.
 * Ensuring such things are at the user's responsability.
 */
class PhyFunct
{
  public:
/*!\brief Default constructor*/
    PhyFunct(){}
/*!\brief Copy constructor*/
    PhyFunct(const PhyFunct &rhs){*this = rhs;}
/*!\brief Default destructor*/
    ~PhyFunct(){}
/*\brief Abscissa*/
    ParameterPhy x;
/*\brief Ordinate*/
    ParameterPhy y;
/*\brief Name*/
    std::string name;

/*\brief showAll()*/
    void showAll(std::ostream &out = std::cout) const;

/*\brief Operator =*/
    PhyFunct &operator=(const PhyFunct &rhs)
                        {
                           x = rhs.x;
                           y = rhs.y;
                           name = rhs.name;
                           return *this;
                        }
};

}

#endif
