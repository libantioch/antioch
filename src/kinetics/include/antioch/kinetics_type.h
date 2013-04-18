//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#ifndef _ANTIOCH_KINETICS_TYPE_H
#define _ANTIOCH_KINETICS_TYPE_H

namespace Antioch{
/*!
 *
 * \class KineticsType
 * \brief base class for kinetics models
 */
class KineticsType{
   public:
      KineticsType(){}
      virtual ~KineticsType(){}

    //! \return the rate evaluated at \p T.
    template <typename StateType>
    virtual StateType operator()(const StateType& T) const = 0.;

    //! \return the derivative with respect to temperature evaluated at \p T.
    template <typename StateType>
    virtual StateType derivative( const StateType& T ) const = 0.;

    //! Simultaneously evaluate the rate and its derivative at \p T.
    template <typename StateType>
    virtual void rate_and_derivative(const StateType& T, StateType& rate, StateType& drate_dT) const = 0.;
}

}

#endif
