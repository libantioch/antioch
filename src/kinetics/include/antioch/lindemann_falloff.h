//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#ifndef _ANTIOCH_LINDEMANN_FALLOFF_H
#define _ANTIOCH_LINDEMANN_FALLOFF_H

//Antioch

//C++

namespace Antioch{

template <typename CoeffType>
class LindemannFalloff{
     public:
       LindemannFalloff();
       ~LindemannFalloff();

     CoeffType operator()(const CoeffType& T) const;

};
  template<typename CoeffType>
  inline
  unsigned int Reaction<CoeffType>::operator(const CoeffType &T) const
  {
    return 1.;
  }


}
