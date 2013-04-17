//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#ifndef _CORE_UNCERTAINTY_
#define _CORE_UNCERTAINTY_

//Antioch

//C++
#include <vector>
#include <string>
#include <cmath>

namespace Antioch{
/****************
 * Definitions  *
 ****************/

#define CORE_UNCERTAINTY_TYPE_ABSOLUTE "absolute"
#define CORE_UNCERTAINTY_TYPE_RELATIVE "relative"
#define CORE_UNCERTAINTY_TYPE_UNKNOWN "unknown"
#define CORE_UNCERTAINTY_TYPE_NONE "none"

/*! \class CoreUnc
 * \brief Class CoreUnc is for uncertainty dealing, at the simplest level
 *
 * Minimal uncertainty is a std::vector of double for values and a std::string for  
 * the type:                           
 *    - unknown (std::string CORE_UNCERTAINTY_TYPE_UNKNOWN)
 *    - none (std::string CORE_UNCERTAINTY_TYPE_NONE)
 *    - absolute (std::string CORE_UNCERTAINTY_TYPE_ABSOLUTE)
 *    - relative (std::string CORE_UNCERTAINTY_TYPE_RELATIVE)
 *
 * This class deals with parameters that have uncertainty of the
 * type \f$\nu \pm \sigma\f$, with \f$\sigma\f$ being the standard
 * uncertainty.
 *
 * Uncertainty calculations are performed following the GUM recommendations.
 * \f[
 *    u_c^2 = \sum_{x_i} \left(\frac{\partial f}{\partial x_i}\right)^2 u_{x_i}^2
 * \f]
 * Thus the class itself can only combined uncertainties for additions and substractions
 * following the previous equation.
 *
 * This class stores values, which are interpreted as the variance
 * (\f$\sigma^2\f$) in case the type is CORE_UNCERTAINTY_TYPE_ABSOLUTE.
 * In this sense, the addition is possible only if the uncertainties are
 * absolute (substraction being the very same operation). The methods
 * containing Variance or StdUnc make the necessary adjustements for
 * each uncertainty type. The methods containing Dvalue are simply
 * dealing with the value and are thus faster, but the user should always
 * be aware of how the value is interpreted when using getters and setters:
 *   - variance in case of type CORE_UNCERTAINTY_TYPE_ABSOLUTE,
 *   - standard uncertainty in case of type CORE_UNCERTAINTY_TYPE_UNKNOWN,
 *   - relative standard uncertainty in case of type CORE_UNCERTAINTY_TYPE_RELATIVE,
 *   - it does not matter in case of type CORE_UNCERTAINTY_TYPE_NONE (as it is zero).
 *
 * The uncertainty type is set to absolute if not provided, this is (and should
 * always be) the default mode.
 */
class CoreUnc{
  public:
/*! \brief Copy constructor, using the operator=*/
   CoreUnc(const CoreUnc &ToBeCopied){*this = ToBeCopied;}
/*! \brief Building constructor, using a std::vector<double> of values and a std::string for the type*/
   CoreUnc(std::vector<double> dval,std::string typeStr = CORE_UNCERTAINTY_TYPE_ABSOLUTE):dvalues(dval),type(typeStr){}
/*! \brief Building constructor, using a double value and a std::string for the type*/
   CoreUnc(double dval, std::string typeStr = CORE_UNCERTAINTY_TYPE_ABSOLUTE);
/*! \brief Building constructor, overload of the previous, only to cover possibilities. No default values in this one*/
   CoreUnc(std::string typeStr, double dval):type(typeStr){addDvalue(dval);}
/*! \brief Building constructor, storing only a type, no values, unknown by default*/
   CoreUnc(std::string typeStr = CORE_UNCERTAINTY_TYPE_UNKNOWN):type(typeStr){}
/*! \brief Default destructor*/
   virtual ~CoreUnc(){}

/*! \brief Values setter.*/
   void setDvalues(const std::vector<double> &dval)    {dvalues = dval;}
/*! \brief Values adder. A std::vector of values is added to the already stored values*/
   void addDvalues(const std::vector<double> &dval);
/*! \brief Values getter. Provide a std::vector<double> of the stored values*/
   const std::vector<double> getDvalues()        const {return dvalues;}
/*! \brief Values getter. Provide a std::vector<double> of the square root of the stored values*/
   const std::vector<double> getSqrtvalues() const;
/*! \brief Set all values to zeros.*/
   void setToNullUnc();

/*! \brief Value getter*/
   double getDvalue(int i)                const;
/*! \brief Overload of double getDvalue(int)*/
   double getDvalue(unsigned int i)       const  {return getDvalue((int)i);}
/*! \brief Value getter, but for the variances*/
   double getSqrtvalue(int i)             const;
/*! \brief Overload of double getSqrtvalue(int)*/
   double getSqrtvalue(unsigned int i)    const {return getSqrtvalue((int)i);}
/*! \brief Value setter*/
   void setDvalue(double dval, int i = 0);
/*! \brief Overload of void setDvalue(double,int)*/
   void setDvalue(double dval, unsigned int i)  {setDvalue(dval,(int)i);}
/*! \brief Value adder. */
   void addDvalue(double dval)                  {dvalues.push_back(dval);}
/*! \brief Square root of value adder. */
   void addSqrtvalue(double dval)               {dvalues.push_back(std::sqrt(dval));}

/*! \brief Size getter. Returns the numbers of stored values*/
   const int sizeDVal()                   const {return dvalues.size();}

/*! \brief Type setter*/
   void setType(const std::string &typeStr)     {type = typeStr;}
/*! \brief Type getter*/
   const std::string getType()            const {return type;}
/*! \brief  Set the type to unknown and clear the values
 *
 * If for any reason we lose the type of uncertainty.
 * Losing an uncertainty can be done with a combination with 
 * an unknown uncertainty, the type can be set to unknown.
 */
   void setTypeToUnknown();
/*! \brief Set the type to none and clear the values
 */
   void setTypeToNone();

/*! \brief Clear the values and the type.*/
   void clear(){dvalues.clear();type.clear();}
/*! \brief Values resizer*/
   void resizeDValues(unsigned int newSize, double dv = 0.) {dvalues.resize(newSize,dv);}
/*! \brief Overload of void resizeValues(unsigned int)*/
   void resizeDValues(int newSize, double dv = 0.)          {resizeDValues((unsigned int)newSize,dv);}
/*! \brief Values reserver*/
   void reserveDValues(unsigned int newSize)                {dvalues.reserve(newSize);}
/*! \brief Overload of void reserveValues(unsigned int)*/
   void reserveDValues(int newSize)                         {reserveDValues((unsigned int)newSize);}


/*! \brief Uncerainties are equal if and only if the type and all the values are
 * the same, in the same order*/
   bool const operator==(const CoreUnc &rhs) const;
/*! \brief Not bool const operator==(const CoreUnc&)*/
   bool const operator!=(const CoreUnc &rhs) const {return (!(*this == rhs));}
/*! \brief Assignement operator, equalizing type and values.*/
   CoreUnc &operator =(const CoreUnc &rhs);
/*! \brief Adding operator
 *
 * At the CoreUnc level, we combine in this operator only
 * uncertainties of type CORE_UNCERTAINTY_TYPE_ABSOLUTE and
 * CORE_UNCERTAINTY_TYPE_NONE. The checks performed are first
 * on the types of the uncertainties, then a second test on
 * the number of values is performed: if the uncertainties are
 * both of type CORE_UNCERTAINTY_TYPE_ABSOLUTE but with different 
 * numbers of uncertainty values, it produces and error.
 *
 * Then the uncertainties are added following
 * the GUM recommandations:
 * \f[
 *   u^2_{add} = u_1^2 + u_2^2
 * \f]
 * */
   CoreUnc &operator+=(const CoreUnc &rhs);
/*! \brief Substracting operator
 *
 * Following the GUM recommandation and as the variance
 * is a squared value:
 * \f[
 *    u_{sub}^2 = u_{add}^2
 * \f]
 * See CoreUnc & operator+=(const CoreUnc&) for more details.
 * */
   CoreUnc &operator-=(const CoreUnc &rhs);
/*! \brief Multiplying operator with a double
 *
 * The operator first check if the type is unknown 
 * (CORE_UNCERTAINTY_TYPE_UNKNOWN or empty std::string).
 * If this is the case, it sends back an error. If not,
 * the calcuations are performed following the GUM 
 * recommandations for uncertainties combinations, which
 * means that the values are multiplied by the square root
 * of the double, if the uncertainty type is CORE_UNCERTAINTY_TYPE_ABSOLUTE
 * or the value if the uncertainty type is CORE_UNCERTAINTY_TYPE_RELATIVE.
 */
   CoreUnc &operator*=(double r);
/*! \brief Dividing operator with a double
 *
 * This is multiplying by the inverse of the double. See
 * CoreUnc &operator*=(double).
 * */
   CoreUnc &operator/=(double r);

/*! \brief Alternative call to CoreUnc &operator=(CoreUnc&)*/
   void equalize(const CoreUnc &rhs)   {*this  = rhs;}
/*! \brief Alternative call to CoreUnc &operator+=(CoreUnc&)*/
   void add(const CoreUnc &rhs)        {*this += rhs;}
/*! \brief Alternative call to CoreUnc &operator-=(CoreUnc&)*/
   void substract(const CoreUnc &rhs)  {*this -= rhs;}
/*! \brief Alternative call to CoreUnc &operator*=(double)*/
   void multiply(double r)             {*this *= r;}
/*! \brief Alternative call to CoreUnc &operator/=(double)*/
   void divide(double r)               {*this /= r;}

  protected:
/*!\brief A std::vector to store the values*/
   std::vector<double> dvalues;
/*!\brief A std::string to store the type*/
   std::string type;
};
}
#endif
