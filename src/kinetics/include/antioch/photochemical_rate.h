//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef _PHOTOCHEMICAL_RATE_
#define _PHOTOCHEMICAL_RATE_

namespace Antioch{

template<typename StateType, typename VectorStateType>
class PhotoRate{
     public:
       PhotoRate(const VectorStateType &cs, const VectorStateType &lambda):
                _cross-section(cs),_lambda_grid(lambda){}
       ~PhotoRate(){}
       StateType forward_rate_constant(const VectorStateType &hv, const VectorStateType &lambda);

       void set_cross_section(const VectorStateType &cs) {_cross-section = cs;}
       void set_lambda_grid(const VectorStateType &l)    {_lambda_grid = l;}

     private:
       VectorStateType _cross-section;
       VectorStateType _lambda_grid;

       VectorStateType make_grid_cool(const VectorStateType &hv, const VectorStateType &lambda);

};

template<typename StateType, typename VectorStateType>
StateType PhotoRate<VectorStateType,StateType>::forward_rate_constant(const VectorStateType &hv, const VectorStateType &lambda)
{
// lambda grid
  VectorStateType hvgrid = this->make_grid_cool(hv,lambda);
// integration, those are bins => just multiply
  StateType rfwd;
  Antioch::set_zero(rfwd);
  for(unsigned int i = 0; i < hvgrid.size(); i++)
  {
    rfwd += _cross_section[i] * hvgrid[i];
  }

  return rfwd;
}

template <typename StateType, typename VectorStateType>
VectorStateType PhotoRate::make_grid_cool(const VectorStateType &hv, const VectorStateType &lambda)
{
  VectorStateType out_hv_on_grid;
  unsigned int j(0);

  for(unsigned int i = 0; i < _lambda_grid.size()-1; i++)//bin per bin, bin hv[j] between lambda[j] and lambda[j+1]
  {
     while(lambda[j] < _lamdba_grid[i])
     {
       j++;
       if(j >= lambda.size())return out_hv_on_grid;
     }
     StateType bin;
     Antioch::set_zero(bin);
     StateType diff_min = lambda[j] - _lambda_grid[i];
     StateType diff_max = _lambda_grid[i+1] - lambda[j];
     bin += hv[j]   * (diff_min)/(lambda[j] - lambda[j-1]);
     if(j < lambda.size() - 1)bin += hv[j+1] * (diff_max)/(lambda[j+1] - lambda[j]);
     out_hv_on_grid.push_back(bin);
  }

  return out_hv_on_grid;
}

}

#endif
