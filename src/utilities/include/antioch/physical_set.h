//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_PHYSICAL_SET_H
#define ANTIOCH_PHYSICAL_SET_H

namespace
{
  template<typename Physics, typename Mixture>
  PhysicalSet
  {
      public:
        PhysicalSet(const Mixture & mix):_mixture(mix)
          {
            physical_set_initialize(_set, physical_set_tag<Physics>::type());
          }

          ~PhysicalSet()
          {
            physical_set_delete(set, physical_set_tag<Physics>::type());
          }

        const Mixture & mixture() const {return _mixture;}

        template <typename InitType>
        void add_model(const std::string & species, const InitType & initMe)
        {
           
           antioch_assert( _mixture.species_name_map().find(species_name) !=
	        	   _mixture.species_name_map().end() );

           unsigned int s = _mixture.species_name_map().find(species)->second;

           physical_set_add(_set,initMe,physical_set_tag<Physics>::type());
        }

        template <typename InitType>
        void add_model(const InitType & initMe)
        {
           physical_set_add(_set,initMe,physical_set_tag<Physics>::type());
        }

        template <typename InitType>
        void reset_model(unsigned int s, const InitType & coefs)
        {
           physical_set_reset(s,_set,initMe,physical_set_tag<Physics>::type());
        }

        template <typename InitType>
        void reset_model(const InitType & coefs)
        {
           physical_set_reset(_set,initMe,physical_set_tag<Physics>::type());
        }

        template<typename StateType>
        ANTIOCH_AUTO(StateType) operator()(unsigned int s, const StateType & T) const 
        ANTIOCH_AUTOFUNC(StateType, physical_set_operator_viscosity(_set,s,T,physical_tag<Physics>::type()))

        template<typename StateType, typename MatrixStateType>
        void operator()(const StateType & T, const StateType & cTot, MatrixStateType & Ds) const 
        {
           physical_set_operator_diffusion_comes_first(_set, T, cTot, Ds, physical_tag<Physics>::type());
        }

        template<typename StateType>
        void operator()(unsigned int s, const StateType & mu, const StateType & dss, const StateType & T, const StateType & rho, StateType & k) const
        {
           physical_set_operator_thermal_conductivity(_set, s, mu, dss, T, rho, k, physical_tag<Physics>::type());
        }

        template<typename StateType, typename ThermoEvaluator>
        void operator()(const StateType & cp, const StateType & k, StateType & ds) const
        {
           physical_set_operator_diffusion_comes_last(_set, cp, k, ds, physical_tag<Physics>::type());
        }

      private:
        const Mixture & _mixture;
        SetOrEquation<Physics, is_physical_set<Physics>::value>::type _set;
  };

}

#endif
