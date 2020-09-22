#ifndef INCLUDED_N_MITRAL_KOPELL2009_HPP
#define INCLUDED_N_MITRAL_KOPELL2009_HPP

#include "insilico/core/engine.hpp"
#include "inputcurrent_twomitral.hpp"
#include "noise.hpp"
#include "currents/ILeak_MC.hpp"
#include "currents/INA_MC.hpp"
#include "currents/INAP_MC.hpp"
#include "currents/IKA_MC.hpp"
#include "currents/IKS_MC.hpp"
#include "currents/IKfast_MC.hpp"



namespace insilico {

  class N_Mitral_KOPELL2009: public Neuron {
  public:
    void ode_set(state_type &variables, state_type &dxdt, const double t,
		 unsigned index) {
      int v_soma_index = engine::neuron_index(index, "v_soma");
      double v_soma = variables[v_soma_index];
      
      // ODE set
	    
      

      ILeak_MC::current(variables, dxdt, t, index);
      INA_MC::current(variables, dxdt, t, index);
      INAP_MC::current(variables, dxdt, t, index);
      IKA_MC::current(variables, dxdt, t, index);
      IKS_MC::current(variables, dxdt, t, index);
      IKfast_MC::current(variables, dxdt, t, index);
     
     
      double INA_MC = engine::neuron_value(index, "INA_MC");
      double INAP_MC = engine::neuron_value(index, "INAP_MC");
      double IKA_MC = engine::neuron_value(index, "IKA_MC");
      double IKS_MC = engine::neuron_value(index, "IKS_MC");
      double IKfast_MC = engine::neuron_value(index, "IKfast_MC");
      double ILeak_MC = engine::neuron_value(index, "ILeak_MC");	
      
    	double noisefactor=engine::neuron_value(index,"noise");
	    double dc=engine::neuron_value(index,"dc");	
	    pulse::inject(index,t);
      noise::inject(index,t);
      double current = engine::neuron_value(index,"currval");
      double variation=engine::neuron_value(index,"noiseval");
      double kacurrent=engine::neuron_value(index,"kacurrent");
      double variate=variation*current*noisefactor;
      double inh_scale=engine::neuron_value(index,"scale");

     
    
      // Input
      //Synapse associated delay for next spikes
      std::vector<unsigned> s_indices = engine::get_pre_neuron_indices(index, "s");
      std::vector<double> esyn_values = engine::get_pre_neuron_values(index, "esyn");
      std::vector<double> gsyn_values = engine::get_pre_neuron_values(index,"g");
      double I_Syn_1=0;

      for(std::vector<int>::size_type iterator = 0; iterator < s_indices.size(); ++iterator) {
	    I_Syn_1 = I_Syn_1 + gsyn_values[iterator]*variables[s_indices[iterator]] *(v_soma - esyn_values[iterator]);
		
	    	    }
      
      //Input current as a synapse
		
      engine::neuron_value(index,"curr",current-variate);
      //engine::neuron_value(index,"Isyn",I_Syn_1*inh_scale);
      engine::neuron_value(index,"variate",variation);
      engine::neuron_value(index,"INA",INA_MC);
      engine::neuron_value(index,"INAP",INAP_MC);
      engine::neuron_value(index,"IKfast",IKfast_MC);
      engine::neuron_value(index,"IKA",IKA_MC);
      engine::neuron_value(index,"IKS",IKS_MC);
      

      	 

      //Injector
	
      dxdt[v_soma_index]=-ILeak_MC-INA_MC-INAP_MC-IKA_MC*kacurrent-IKS_MC-IKfast_MC+current-variate-I_Syn_1;//-I_Syn_1;//+current-variate-I_Syn_1;
    } // function ode_set
  }; // class N_SquidAxon_HH1952

} // insilico
#endif


