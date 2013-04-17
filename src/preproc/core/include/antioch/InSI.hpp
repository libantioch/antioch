//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
#ifndef _IN_SI_
#define _IN_SI_

//Antioch

//C++
#include <iostream>

namespace Antioch{

/*! \class InSI
 * \brief Seven integers to characterize the power vector
 *
 * A unit is caracterize by the power associated to
 * every base unit, this class stores those power.
 */
class InSI{
    public:
/*! \brief Copy constructor, uses InSI &operator=(const InSI&)*/
      InSI(InSI const &rhs){*this = rhs;}
/*! \brief Building constructor, fully descriptive, with zeros as default values to
 * simplify coder's interface*/
      InSI(int i0=0,int i1=0, int i2=0, int i3=0, int i4=0, int i5=0, int i6=0, int i7=0):
        m(i0),kg(i1),s(i2),A(i3),K(i4),mol(i5),cd(i6),rad(i7){}

/*! \brief << operator, to format the power vector*/
      friend std::ostream &operator<< (std::ostream &out, const InSI & rhs)
        {
          out << "["
              << rhs.get_m()   << "(m),"
              << rhs.get_kg()  << "(kg),"
              << rhs.get_s()   << "(s),"
              << rhs.get_A()   << "(A),"
              << rhs.get_K()   << "(K),"
              << rhs.get_mol() << "(mol),"
              << rhs.get_cd()  << "(cd),"
              << rhs.get_rad() << "(rad)]";
          return out;
        }

/*! \brief Bool equalize operator, true if all the powers are equal*/
      bool const operator== (const InSI & rhs) const;
/*! \brief Not bool const operator==(const InSI&)*/
      bool const operator!= (const InSI & rhs) const;
/*! \brief Assignement operator, equalize all the powers*/
      InSI &     operator=  (const InSI & rhs);
/*! \brief Adding operator, add all the powers*/
      InSI &     operator+= (const InSI & rhs);
/*! \brief Substracting operator, substract all the powers*/
      InSI &     operator-= (const InSI & rhs);
/*! \brief Multiplying operator, multiply all the powers*/
      InSI &     operator*= (int rhs);
/*! \brief Dividing operator.
 *
 * Dividing a power needs first to check that the
 * division is possible, i.e. if the power is not
 * null, it must be a multiple of the divider. If not
 * it sends back an error and stops. This division
 * exist for root purposes.
 */
      InSI &     operator/= (int rhs);
/*! \brief Adding operator, add all the powers*/
      InSI       operator+  (const InSI & rhs) const;
/*! \brief Substracting operator, substract all the powers*/
      InSI       operator-  (const InSI & rhs) const;
/*! \brief Multiplying operator, multiply all the powers*/
      InSI       operator*  (int rhs) const;
/*! \brief Dividing operator, see InSi & operator/=(int) for details*/
      InSI       operator/  (int rhs) const;

/*! \brief Set all the powers to zeros*/
      void clear();

/*! \brief meter power getter*/
      const int get_m()   const {return m;}
/*! \brief kilogramme power getter*/
      const int get_kg()  const {return kg;}
/*! \brief second power getter*/
      const int get_s()   const {return s;}
/*! \brief ampere power getter*/
      const int get_A()   const {return A;}
/*! \brief kelvin power getter*/
      const int get_K()   const {return K;}
/*! \brief mol power getter*/
      const int get_mol() const {return mol;}
/*! \brief candela power getter*/
      const int get_cd()  const {return cd;}
/*! \brief radian power getter*/
      const int get_rad() const {return rad;}

/*!\brief Check if empty (all values to zero)*/
      bool empty() const;

    private:
      int m,kg,s,A,K,mol,cd,rad;
};

}

#endif
