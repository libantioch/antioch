//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_HERCOURT-ESSEN_RATE_H
#define ANTIOCH_HERCOURT-ESSEN_RATE_H

// C++
#include <cmath>
#include <iostream>

namespace Antioch
{
  //! Hercourt-Essen rate equation.
  /*!
   * Hercourt-Essen rate equation.  Computes rates of the form
   * \f$ C_f\times T^\eta \f$.
   */
  template<typename CoeffType=double>
  class HercourtEssenRate
  {
  
  public:

    HercourtEssenRate (const CoeffType Cf=0., const CoeffType eta=0., const CoeffType Ea=0.);
    ~HercourtEssenRate();
    
    void set_Cf( const CoeffType Cf );
    void set_eta( const CoeffType eta );

    CoeffType Cf() const;
    CoeffType eta() const;

    //! \return the rate evaluated at \p T.
    template <typename StateType>
    StateType operator()(const StateType& T) const;

    //! \return the derivative with respect to temperature evaluated at \p T.
    template <typename StateType>
    StateType derivative( const StateType& T ) const;

    //! Simultaneously evaluate the rate and its derivative at \p T.
    template <typename StateType>
    void rate_and_derivative(const StateType& T, StateType& rate, StateType& drate_dT) const;

    //! Formatted print, by default to \p std::cout
    void print(std::ostream& os = std::cout) const;

    //! Formatted print.
    friend std::ostream& operator<<(std::ostream& os, const HercourtEssenRate& rate)
    {
      rate.print(os);
      return os;
    }

  private:

    CoeffType _Cf;
    CoeffType _eta;
    
  };

  template<typename CoeffType>
  HercourtEssenRate<CoeffType>::HercourtEssenRate(const CoeffType Cf, const CoeffType eta)
    : _Cf(Cf),
      _eta(eta),
  {
    return;
  }

  template<typename CoeffType>
  HercourtEssenRate<CoeffType>::~HercourtEssenRate()
  {
    return;
  }

  template<typename CoeffType>
  void HercourtEssenRate<CoeffType>::print(std::ostream& os) const
  {
    os << _Cf;
    os << "*T^" << _eta;

    return;
  }

  /* ------------------------- Inline Functions -------------------------*/
  template<typename CoeffType>
  inline
  void HercourtEssenRate<CoeffType>::set_Cf( const CoeffType Cf )
  {
    _Cf = Cf;
    return;
  }

  template<typename CoeffType>
  inline
  void HercourtEssenRate<CoeffType>::set_eta( const CoeffType eta )
  {
    _eta = eta;
    return;
  }

  template<typename CoeffType>
  inline
  CoeffType HercourtEssenRate<CoeffType>::Cf() const
  { return _Cf; }

  template<typename CoeffType>
  inline
  CoeffType HercourtEssenRate<CoeffType>::eta() const
  { return _eta; }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType HercourtEssenRate<CoeffType>::operator()(const StateType& T) const
  {
    using std::pow;
    return _Cf* (pow(T,_eta));
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  StateType HercourtEssenRate<CoeffType>::derivative( const StateType& T ) const
  {
    return (*this)(T)/T*(_eta);
  }

  template<typename CoeffType>
  template<typename StateType>
  inline
  void HercourtEssenRate<CoeffType>::rate_and_derivative( const StateType& T,
						      StateType& rate,
						      StateType& drate_dT) const
  {
    rate     = (*this)(T);
    drate_dT = rate/T*(_eta);
    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_ARRHENIUS_RATE_H
