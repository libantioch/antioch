//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_SIGMA_BIN_MANAGER_H
#define ANTIOCH_SIGMA_BIN_MANAGER_H

//Antioch
#include "antioch/vector_utils.h"

//Planet

//C++
#include <vector>
#include <iostream>

namespace Antioch{

template <typename VectorCoeffType = std::vector<double> >
class SigmaBinConverter
{
     public:
        SigmaBinConverter();
        ~SigmaBinConverter();

        void y_on_custom_grid(const VectorCoeffType &x_old, const VectorCoeffType &y_old,  const VectorCoeffType &x_new);

     private:
};

template <typename VectorCoeffType>
inline
SigmaBinConverter<VectorCoeffType>::SigmaBinConverter()
{
  return;
}

template <typename VectorCoeffType>
inline
SigmaBinConverter<VectorCoeffType>::~SigmaBinConverter()
{
  return;
}

template <typename VectorCoeffType>
inline
void SigmaBinConverter<VectorCoeffType>::y_on_custom_grid(const VectorCoeffType &x_old, const VectorCoeffType &y_old,  
                                                          const VectorCoeffType &x_custom,    VectorCoeffType &y_custom)
{
  
  y_custom.clear();
  y_custom.resize(x_custom.size(),0.L);

  unsigned int ilow(0);
  while(x_custom[ilow] <= x_old[0])ilow++; //skipping too low bins

  unsigned int j(1);
  for(unsigned int i = ilow; i < x_custom.size() - 1; i++)//bin per bin, right stairs: y_old[i] from x_old[i] to x_old[i+1]
  {
     while(x_old[j] < x[i]) //find lowest j / x_old[j] > x[i]
     {
       j++;
      if(j >= x_old.size() - 1)return;
     }
     typedef typename value_type<VectorCoeffType>::type CoeffType bin;
     Antioch::set_zero(bin);

//here we are: x_old[j-1] =< x_custom[i] < x_old[j] with j-1 >= 0

     //targeted bin within stored bin: x_old[j-1] =< x_custom[i] < x_custom[i+1] =< x_old[j], take all of them
     while(x_old[j] >= x_custom[i+1])
     {
       y_custom[i] = y_old[j-1]; //rectangle from i to i+1, same height
       i = i + 1;
       if(i >= x_custom.size() - 2)break;
     }

     //x_old[j-1] < x_custom[i] < x_old[j] < x_custom[i+1], calculating rectangle from x_custom[i] to x_old[j], height is y_old[j-1]
     bin = y_old[j-1] * (x_old[j] - x_custom[i]); 

// finding lowest j / x_custom[i+1] < x_old[j], 
// adding all the k cases x_old[j-1] < x_custom[i] < x_old[j] < x_old[j+1] < ... < x_old[j+k] < x_custom[i+1]
     j++;
     while(x_old[j] < x_custom[i+1])
     {
        bin += y_old[j-1] * (x_old[j] - x_old[j-1]); // adding contained bins
        j++;
        if(j >= x_old.size()-1)break;
     }

//now we have found k_max, we calculate the rectangle from x_old[j + k_max] to x_custom[i+1]
     if(j < x_old.size() - 1)bin += y_old[j-2] * (x_custom[i+1] - x_old[j-1]); //if exist, above rectangle
     y_custom[i] = bin/(x_custom[i+1] - x_custom[i]); //rectangle from i to i+1
  }

  return;
}

}

#endif
