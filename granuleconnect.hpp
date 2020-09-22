#ifndef INCLUDED_GC_HPP
#define INCLUDED_GC_HPP

#include "insilico/core.hpp"

#include "IKA_GC.hpp"
#include "IKM_GC.hpp"
#include "INA_GC.hpp"
#include "IK_GC.hpp"
#include "ILeak_GC.hpp"



namespace insilico {

class gc: public Neuron {
 public:
   void ode_set(state_type &variables, state_type &dxdt, const double t,
               unsigned index) {
    int v_soma_index = engine::neuron_index(index, "v_soma");
    double v_soma = variables[v_soma_index];
    int v_dend_index = engine::neuron_index(index, "v_dend");
    double v_dend = variables[v_dend_index];
		
   
        // ODE set
    IKA_GC_soma::current(variables, dxdt, t, index);
    IKM_GC_soma::current(variables, dxdt, t, index);
		INA_GC_soma::current(variables, dxdt, t, index);
    ILeak_GC_soma::current(variables, dxdt, t, index);
    IK_GC_soma::current(variables, dxdt, t, index);
    
    INA_GC_dend::current(variables, dxdt, t, index);
		IK_GC_dend::current(variables, dxdt, t, index);
	  ILeak_GC_dend::current(variables, dxdt, t, index);
			
	  
		double ILeak_soma = engine::neuron_value(index, "ILeak_GC_soma");
    double INA_soma = engine::neuron_value(index, "INA_GC_soma");
    double IK_soma = engine::neuron_value(index, "IK_GC_soma");
  	double IKA_soma = engine::neuron_value(index, "IKA_GC_soma");
  	double IKM_soma = engine::neuron_value(index, "IKM_GC_soma");

	 
	  //double gsyn=engine::neuron_value(index,"g");
	  double IK_dend = engine::neuron_value(index, "IK_GC_dend");
  	double INA_dend = engine::neuron_value(index, "INA_GC_dend");
	  double ILeak_dend = engine::neuron_value(index, "ILeak_GC_dend");
    double dc=engine::neuron_value(index,"dc"); 

	// note the spike
   /* double last_spiked = engine::neuron_value(index, "last_spike");
    double spike_duration = engine::neuron_value(index, "spike_duration");
    double thresh = engine::neuron_value(index, "thresh");

    // associated delay for next spikes
    if((v_dend > thresh) && (t - last_spiked) > spike_duration){
      engine::neuron_value(index, "last_spike", t);
    }
	
 		*/
	  std::vector<unsigned> s_indices = engine::get_pre_neuron_indices(index, "s");
    std::vector<double> esyn_values = engine::get_pre_neuron_values(index, "esyn");
    std::vector<double> gsyn_values = engine::get_pre_neuron_values(index,"g");

		double I_Syn_gc = 0;
		for(std::vector<int>::size_type iterator = 0; iterator < s_indices.size(); ++iterator) {
				I_Syn_gc = I_Syn_gc + gsyn_values[iterator]*variables[s_indices[iterator]] *(v_dend - esyn_values[iterator]);
			
	}
  
	/*engine::neuron_value(index,"Isyn",I_Syn_gc);*/


  //std::cout<<t<<","<<I_Syn_gc<<"\n";
  
  	dxdt[v_dend_index]=-ILeak_dend-IK_dend-INA_dend+(0.01*(v_soma-v_dend))-I_Syn_gc;//-0.0282;//+1
	dxdt[v_soma_index]=-ILeak_soma-IK_soma-INA_soma-IKA_soma-IKM_soma+(0.01*(v_dend-v_soma));//-(1.3866*(v_priden - v_soma)/4);
	/*
	dxdt[v_dend_index]=-ILeak_dend-IK_dend-INA_dend;//-I_Syn_gc;//+(0.01*(v_soma-v_dend));//(1000.241*I_Syn_gc)+(.6822*(v_priden-v_dend)/4);
	dxdt[v_soma_index]=-ILeak_soma-IK_soma-INA_soma-IKA_soma-IKM_soma;//+(0.000001*(v_dend-v_soma));//-(1.3866*(v_priden - v_soma)/4);
*/
  } // function ode_set
}; // class N_SquidAxon_HH19521

}
#endif



