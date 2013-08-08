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
#ifndef _PHYSICAL_ALGEBRA_
#define _PHYSICAL_ALGEBRA_

#include "PhyFunct.hpp"
#include <vector>

namespace Antioch{
/*!
 * \class PhyHarmFunct
 *
 * \brief Code for several sets of \f$y\f$ values for the same set of
 * \f$x\f$ values.
 *
 *
 * This class will take a physical parameter as abscissa and take it as
 * a grid for functions to be given. The given functions must have a
 * \f$x\f$ parameter homogeneous to the abscissa. The \f$y\f$ values
 * will be put on the grid defined by the abscissa, using either zeros
 * values if the gris is outside the \f$x\f$ range of the function, or
 * a simple linear interpolation between the two closest framing values.
 */

class PhyHarmFunct
{
  public:
/*!\brief Default constuctor*/
    PhyHarmFunct(){}
/*!\brief Constuctor given only functions
 *
 * This class provides two ways of defining the common abscissa grid combined
 * with two ways of calculating the value.
 * Defining an abscissa can be done:
 *   - by explicitely giving it,
 *   - by calculation of a regular grid;
 * calculating a value:
 *   - by interpolating between the two surrounding values,
 *   - by calculating the mean value on an interval.
 *
 * By default the definition is the regular grid on a hundred points and
 * the calculation is the interpolation. To chose between different definitions, the 
 * user can provide a key word. The word ``abscissa'' will provoke the
 * constructor to chose a given abscissa, at index xref. Any other word (including no
 * word at all) will cause the constructor to build the regular grid.
 *
 * The calculation method is chosen by the keyword defcal, "interpolation" or
 * "mean".
 *
 * \param functions  set of functions to harmonize
 * \param defcal key to choose between interpolation or mean calculation method
 * \param xref  index of explicit abscissa grid or number or values on
 *                        the regular grid
 * \param deftype  key used to choose between the explicit and regular grid methods.
 *              Default is regular grid, set it to ``abscissa'' to switch to
 *              explicit method.
 */
    PhyHarmFunct(const std::vector<PhyFunct> &functions, const std::string defcal = "interpolation", 
                                                         int xref = -1, const std::string deftype = std::string());
/*!\brief Constructor with xstep
 *
 * The abscissa in this method is given by the min and max
 * found in the functions, and the step xstep.
 */
    PhyHarmFunct(const std::vector<PhyFunct> &functions, double xstep);
/*!\brief Constuctor to set directly an explicit grid*/
    PhyHarmFunct(const ParameterPhy &xgrid, const std::vector<PhyFunct> &functions, const std::string defcal = "interpolation");
/*!\brief Default destructor*/
    ~PhyHarmFunct(){}

/*!\brief Calculate the ordinate on the grid
 *
 * The ordinate of the given function is calculated on
 * the abscissa grid, which must be define before using this
 * method, and the corresponding ParameterPhy is produced. 
 *
 * The defcal string is the keyword between the two calculations methods.
 */
    ParameterPhy setOrdToGrid(const PhyFunct &f, const std::string defcal = "interpolation");

/*!\brief Reference abscissa grid setter*/
    void setAbs(const ParameterPhy &xgrid) {abs = xgrid;}
/*!\brief Add an ordinate, sets the function to the grid*/
    void addOrd(const PhyFunct &f, const std::string defcal = "interpolation");
/*!\brief Add an ordinate from a PhyHarmFunct, sets the ordinate to the grid*/
    void addOrd(const PhyHarmFunct &f, int iord, const std::string defcal = "interpolation");
/*!\brief Add several ordinates*/
    void addOrd(const std::vector<PhyFunct> &f, const std::string defcal = "interpolation");
/*!\brief Add an ordinate, test number of values*/
    void addOrd(const ParameterPhy &f);
/*!\brief Add a zero-filled ordinate with given name and unit*/
    void addZeroOrd(const std::string &name = "Zero ordinate", const std::string &units = std::string());
/*!\brief Test if zero-filled*/
    bool isZeroOrd(int iord) {return hasOnlyZero[iord];}
/*!\brief Test if zero-filled by name*/
    bool isZeroOrd(const std::string &name) {return hasOnlyZero[indexByName(name)];}
/*!\brief Set an ordinate*/
    void setOrd(const PhyFunct &f, int iord = 0, const std::string defcal = "interpolation");
/*!\brief Clear the ordinate vector*/
    void clearOrd()                        {ord.clear();}

/*!\brief Reference abscissa grid getter*/
    ParameterPhy &getAbs()                       {return abs;}
/*!\brief Abscissa grid getter*/
    ParameterPhy  getAbs()                 const {return abs;}
/*!\brief Reference ordinate getter*/
    ParameterPhy &getOrd(int iord = 0)           {return ord[iord];}
/*!\brief Ordinate getter*/
    ParameterPhy  getOrd(int iord = 0)     const {return ord[iord];}
/*!\brief Reference ordinate getter by name*/
    ParameterPhy &getOrd(const std::string &nameOrd) {return ord[indexByName(nameOrd)];}
/*!\brief Ordinate getter by name*/
    ParameterPhy  getOrd(const std::string &nameOrd) const;
/*!\brief Number of ordinates getter*/
    int nOrdinates()                       const {return (int) ord.size();}
/*!\brief Reference vector of ordinates getter*/
    std::vector<ParameterPhy> &getAllOrd()       {return ord;}
/*!\brief Vector of ordinates getter*/
    std::vector<ParameterPhy>  getAllOrd() const {return ord;}

/*!\brief Name getter*/
    std::string getName()                  const {return name;}
/*!\brief Name setter*/
    void setName(const std::string str)          {name = str;}

/*!\brief ShowAll()*/
    void showAll(std::ostream &out = std::cout) const;

/*!\brief Assignement operator.*/
    PhyHarmFunct & operator= (const PhyHarmFunct &rhs);
/*!\brief Egality operator.*/
    bool operator== (const PhyHarmFunct &rhs);

  private:

/*!\brief Index of ordinate of name std::string name*/
    int indexByName(const std::string &nameOrd);

/*!\brief Abscissa calculator
 *
 * This is the method that performs the reference abscissa grid calculations.
 * \param type is the trigger between an arbitrary grid and a regular grid,
 * \param xref is the index of the arbitrary grid or the desired number of
 * values. The default values of the numbers of values on the grid is 100.
 */
    ParameterPhy HarmonizeAbscissa(const std::vector<PhyFunct> *xs, int xref = -1, const std::string type = std::string());
/*!\brief Abscissa calculator
 *
 * This is the method that performs the reference abscissa grid calculations.
 * This version finds the min and max in the vector of PhyFunct, the
 * step is given.
 */
    ParameterPhy HarmonizeAbscissa(const std::vector<PhyFunct> *xs, double xstep);

/*!\brief unit checker
 *
 * It checks if the abscissas are all homogeneous
 */
    bool abscissaHomogeneous(const std::vector<PhyFunct> *xs);

/*!\brief Interpolation method to calculate y*/
    ParameterPhy interpolationOnGrid(const PhyFunct &f);
/*!\brief Integration method to calculate y
 *
 * The integration, mean of values, is calculated between
 * x(i) and x(i+1), last value is extrapolation.
 */
    ParameterPhy integrationOnGrid(const PhyFunct &f);
/*\brief Test if given ParameterPhy is filled with 0.*/
    bool isZeroPP(const ParameterPhy &f);

/*!\brief Interpolation between two data points (x1,y1) (x2,y2)*/
    void interpolation(double &a, double &b, double x1, double y1, double x2, double y2);

/*!\brief A ParameterPhy for the reference abscissa grid*/
    ParameterPhy abs;
/*!\brief A vector of ParameterPhy for the ordinates*/
    std::vector<ParameterPhy> ord;
/*!\brief A string for a name (identifier)*/
    std::string name;
/*!\brief A vector of boolean for zero-filled ordinates to avoid unecessary loops*/
    std::vector<bool> hasOnlyZero;
};
}
#endif
